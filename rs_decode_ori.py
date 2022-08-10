#http://www.xiaojiejing.com/2018/12/07/reed-solomon%e7%a0%81/
gf_exp = [0] * 512 #512个元素的列表. Python 2.6+可以考虑使用bytearray类型
gf_log = [0] * 256

def init_tables(prim=0x11d):
    #使用参数prim给的本原多项式计算指数表和对数表备用
    global gf_exp, gf_log
    gf_exp = [0] * 512 # anti-log (exponential指数) table
    gf_log = [0] * 256 # log(对数) table
    # 计算每一个GF(2^8)域内正整数的指数和对数
    x = 1
    for i in range(0, 255):
        gf_exp[i] = x # 存储指数表
        gf_log[x] = i # 存储对数表
        # 更一般的情况用下面这行，不过速度较慢
        #x = gf_mult_noLUT(x, 2, prim)

        # 只用到 generator==2 或指数底为 2的情况下，用下面的代码速度快过上面的 gf_mult_noLUT():
        x <<= 1 
        if x & 0x100: # 等效于 x >= 256, 但要更快些 (because 0x100 == 256，位运算速度优势)
            x ^= prim # substract the primary polynomial to the current value (instead of 255, so that we get a unique set made of coprime numbers), this is the core of the tables generation

    # Optimization: 双倍指数表大小可以省去为了不出界而取模255的运算 (因为主要用这个表来计算GF域乘法，仅此而已).
    for i in range(255, 512):
        gf_exp[i] = gf_exp[i - 255]
    return [gf_log, gf_exp]

def gf_pow(x, power):
    return gf_exp[(gf_log[x] * power) % 255]

def gf_inverse(x):
    return gf_exp[255 - gf_log[x]] # gf_inverse(x) == gf_div(1, x)

def gf_div(x,y):
    if y==0:
        raise ZeroDivisionError()
    if x==0:
        return 0
    return gf_exp[(gf_log[x] + 255 - gf_log[y]) % 255]

def gf_mul(x,y):
    if x==0 or y==0:
        return 0
    return gf_exp[gf_log[x] + gf_log[y]] # should be gf_exp[(gf_log[x]+gf_log[y])%255] if gf_exp wasn't oversized

def gf_sub(x,y):
    #print(x^y)
    return x^y

def gf_poly_add(p,q):
    r = [0] * max(len(p),len(q))
    for i in range(0,len(p)):
        r[i+len(r)-len(p)] = p[i]
    for i in range(0,len(q)):
        r[i+len(r)-len(q)] ^= q[i]
    return r

def gf_poly_div(dividend, divisor):
    '''适用于GF(2^p)域的快速多项式除法.'''
    # 注意: 多项式系数需要按幂次由高到低排序. 例如: 1 + 2x + 5x^2 = [5, 2, 1], 而非 [1, 2, 5]

    msg_out = list(dividend) # 复制被除数(尾部后缀ecc字节, 用0填充)
    #normalizer = divisor[0] # precomputing for performance
    for i in range(0, len(dividend) - (len(divisor)-1)):
        #msg_out[i] /= normalizer 
        coef = msg_out[i]
        if coef != 0: # 避免log(0)未定义错误.
            for j in range(1, len(divisor)): # 因为多项式的首位都是1, (1^1==0)所以可以跳过
                if divisor[j] != 0: # log(0) is undefined
                    msg_out[i + j] ^= gf_mul(divisor[j], coef) # 等价于数学表达式:msg_out[i + j] += -divisor[j] * coef ,但异或运高效

    # msg_out 包含商和余数, 余数的最高幂次( == length-1)和除数一样, 下面计算分断点.
    separator = -(len(divisor)-1)
    return msg_out[:separator], msg_out[separator:] # 返回商, 余数.

def gf_poly_eval(p,x):
    y = p[0]
    #print("gf_pow",x)
    for i in range(1, len(p)):
        y = gf_mul(y,x) ^ p[i]
    return y

def gf_poly_scale(p,x):
    r = [0] * len(p)
    for i in range(0, len(p)):
        r[i] = gf_mul(p[i], x)
    return r

def rs_calc_syndromes(msg, nsym):
    '''给定原信息和纠错标记数，计算伴随式
    从数学角度看，这个过程就是一个傅里叶变换 (Chien搜索刚好相反).
    '''    
    synd = [0] * nsym
    for i in range(0, nsym):
        #synd[i] = gf_poly_eval(msg, gf_pow(2,i+1)) #XKE
        synd[i] = gf_poly_eval(msg, gf_pow(2,i))
    return [0] + synd # pad with one 0 for mathematical precision (else we can end up with weird calculations sometimes)
    #XKE
    #return synd # pad with one 0 for mathematical precision (else we can end up with weird calculations sometimes)

def rs_generator_poly(nsym):
    g = [1]
    for i in range(0,nsym):
       # g = gf_poly_mul(g, [1, gf_pow(2, i)]) # 指数运算，跟下面等效
       #g = gf_poly_mul(g, [1, gf_exp[i+1]]) #XKE
       g = gf_poly_mul(g, [1, gf_exp[i]])
    return g

def gf_poly_mul(p,q):
    r = [0] * (len(p)+len(q)-1)
    for j in range(0, len(q)):
        for i in range(0, len(p)):
            r[i+j] ^= gf_mul(p[i], q[j])
    return r

def rs_encode_msg(msg_in, nsym):
    '''Reed-Solomon main encoding function, using polynomial division (algorithm Extended Synthetic Division)'''
    if (len(msg_in) + nsym) > 255: raise ValueError("Message is too long (%i when max is 255)" % (len(msg_in)+nsym))
    gen = rs_generator_poly(nsym)
    print("gen poly",gen)
    # 用msg_in值初始化 msg_out 并后缀ecc字节, 用0填充.
    msg_out = [0] * (len(msg_in) + len(gen)-1)
    msg_out[:len(msg_in)] = msg_in

    # Synthetic division main loop
    for i in range(len(msg_in)):
        coef = msg_out[i]

        if coef != 0:  # 避免log(0) 未定义错误
            for j in range(1, len(gen)): # 因为多项式的首位都是1, (1^1==0)所以可以跳过
                msg_out[i+j] ^= gf_mul(gen[j], coef) # equivalent to msg_out[i+j] += gf_mul(gen[j], coef)

    # 除完之后, msg_out 包含商 msg_out[:len(msg_in)]，余数 msg_out[len(msg_in):]. 
    #RS码只用到余数, 所以我们用原信息覆盖商得到完整编码.
    msg_out[:len(msg_in)] = msg_in

    return msg_out

def rs_check(msg, nsym):
    '''Returns true if the message + ecc has no error, false otherwise '''
    return ( max(rs_calc_syndromes(msg, nsym)) == 0 )

def rs_find_error_locator(synd, nsym, erase_loc=None, erase_count=0):
    '''Find error/errata locator and evaluator polynomials with Berlekamp-Massey algorithm'''
    # 使用BM算法迭代计算error locator polynomial.计算差异项 Delta, 由此判断是否更新 error locator polynomial.

    # Init the polynomials
    if erase_loc: 
        err_loc = list(erase_loc)
        old_loc = list(erase_loc)
    else:
        err_loc = [1] # Sigma(errors/errata locator polynomial) .
        old_loc = [1] # BM 是一个迭代算法，需要对比新旧值(Sigma)来判断哪些变量需要更新。

    synd_shift = 0
    if len(synd) > nsym: synd_shift = len(synd) - nsym

    for i in range(0, nsym-erase_count): 
        if erase_loc: # if an erasures locator polynomial was provided to init the errors locator polynomial, then we must skip the FIRST erase_count iterations (not the last iterations, this is very important!)
            K = erase_count+i+synd_shift
        else: # if erasures locator is not provided, then either there's no erasures to account or we use the Forney syndromes, so we don't need to use erase_count nor erase_loc (the erasures have been trimmed out of the Forney syndromes).
            K = i+synd_shift

        # Compute the discrepancy Delta
        # 计算 Delta: error locator和 the syndromes进行多项式乘法运算。
        #delta = gf_poly_mul(err_loc[::-1], synd)[K] # 严格的讲应该是 gf_poly_add(synd[::-1], [1])[::-1] , 但对于正确解码而言并不必要.
        # 而且同时还能优化效率: 我们只计算第K项的多项式乘法. 避免套嵌循环.
        delta = synd[K]
        for j in range(1, len(err_loc)):
            delta ^= gf_mul(err_loc[-(j+1)], synd[K - j]) 
        #print("delta", K, delta, list(gf_poly_mul(err_loc[::-1], synd))) # debugline

        # Shift polynomials to compute the next degree
        old_loc = old_loc + [0]

        # 计算 the errata locator and evaluator polynomials
        if delta != 0: # delta非0才更新
            if len(old_loc) > len(err_loc): 
                # 计算 Sigma(errata locator polynomial)
                new_loc = gf_poly_scale(old_loc, delta)
                old_loc = gf_poly_scale(err_loc, gf_inverse(delta)) # err_loc * 1/delta = err_loc // delta
                err_loc = new_loc

            # 更新Delta
            err_loc = gf_poly_add(err_loc, gf_poly_scale(old_loc, delta))

    # 查看结果是否正确, 或者错误太多无法修复
    while len(err_loc) and err_loc[0] == 0: del err_loc[0] # 排除前面的0
    errs = len(err_loc) - 1
    if (errs-erase_count) * 2 + erase_count > nsym:
        #raise ReedSolomonError("Too many errors to correct")    # 错误太多
        print("Too many errors to correct")    # 错误太多

    return err_loc

def rs_find_errors(err_loc, nmess): # nmess is len(msg_in)
    '''Find the roots (ie, where evaluation = zero) of error polynomial by brute-force trial, this is a sort of Chien's search
    (but less efficient, Chien's search is a way to evaluate the polynomial such that each evaluation only takes constant time).'''
    errs = len(err_loc) - 1
    err_pos = []
    for i in range(nmess): # normally we should try all 2^8 possible values, but here we optimize to just check the interesting symbols
        if gf_poly_eval(err_loc, gf_pow(2, i)) == 0: # It's a 0? Bingo, it's a root of the error locator polynomial,
                                                     # in other terms this is the location of an error
            err_pos.append(nmess - 1 - i)
    # Sanity check: the number of errors/errata positions found should be exactly the same as the length of the errata locator polynomial
    if len(err_pos) != errs:
        # couldn't find error locations
        #raise ReedSolomonError("Too many (or few) errors found by Chien Search for the errata locator polynomial!")
        print("Too many (or few) errors found by Chien Search for the errata locator polynomial!")
    return err_pos

def rs_forney_syndromes(synd, pos, nmess):
    # Compute Forney syndromes, which computes a modified syndromes to compute only errors (erasures are trimmed out). 别混淆Forney syndromes 和 Forney algorithm(根据错误位置纠正信息的算法).
    erase_pos_reversed = [nmess-1-p for p in pos] # prepare the coefficient degree positions (instead of the erasures positions)

    # Optimized method, all operations are inlined
    fsynd = list(synd[1:])      # make a copy and trim the first coefficient which is always 0 by definition
    for i in range(0, len(pos)):
        x = gf_pow(2, erase_pos_reversed[i])
        for j in range(0, len(fsynd) - 1):
            fsynd[j] = gf_mul(fsynd[j], x) ^ fsynd[j + 1]

    # 上述算法是进过优化之后的, 等价的理论算法应该是: fsynd = (erase_loc * synd) % x^(n-k)
    #erase_loc = rs_find_errata_locator(erase_pos_reversed, generator=generator) # computing the erasures locator polynomial
    #fsynd = gf_poly_mul(erase_loc[::-1], synd[1:]) # then multiply with the syndrome to get the untrimmed forney syndrome
    #fsynd = fsynd[len(pos):] # 删除(trim)没用的前 erase_pos 个系数. 似乎没必要, 不过能降低后面BM算法的计算量.

    return fsynd

def rs_find_errata_locator(e_pos):
    '''误码定位
       eg："_"为误码，"h_ll_ worldxxxxxxxxx"中误码位置应为： n-1 - [1, 4] = [18, 15] = erasures_loc.'''

    e_loc = [1] # 初始化为1而非0是因为要计算乘法
    # erasures_loc = product(1 - x*alpha**i) for i in erasures_pos and where alpha is the alpha chosen to evaluate polynomials.
    for i in e_pos:
        e_loc = gf_poly_mul( e_loc, gf_poly_add([1], [gf_pow(2, i), 0]) )
        #XKE
        #e_loc = gf_poly_mul( e_loc, gf_poly_add([1], [gf_pow(2, i+1), 0]) )
    return e_loc

def rs_find_error_evaluator(synd, err_loc, nsym):
    '''Compute the error (or erasures if you supply sigma=erasures locator polynomial, or errata) evaluator polynomial Omega
       from the syndrome and the error/erasures/errata locator Sigma.'''

    # Omega(x) = [ Synd(x) * Error_loc(x) ] mod x^(n-k+1)
    _, remainder = gf_poly_div( gf_poly_mul(synd, err_loc), ([1] + [0]*(nsym+1)) ) # 除法运算只是为了截短

    # Faster way that is equivalent
    #remainder = gf_poly_mul(synd, err_loc) # 乘法
    #remainder = remainder[len(remainder)-(nsym+1):] # 截短

    return remainder

def rs_correct_errata(msg_in, synd, err_pos): # err_pos is a list of the positions of the errors/erasures/errata
    '''Forney algorithm, computes the values (error magnitude) to correct the input message.'''
    # 结合前面的定位结果进行错误判别(errata locator polynomial)
    coef_pos = [len(msg_in) - 1 - p for p in err_pos] # 位置换算为多项式系数 (eg: 由 [0, 1, 2] 换算为 [len(msg)-1, len(msg)-2, len(msg) -3])
    err_loc = rs_find_errata_locator(coef_pos)
    # 计算 errata evaluator polynomial (术语叫 Omega 或 Gamma)
    err_eval = rs_find_error_evaluator(synd[::-1], err_loc, len(err_loc)-1)[::-1]

    # Chien 算法搜索 err_pos 得到误码位置 X  (the roots of the error locator polynomial, ie, where it evaluates to 0)
    X = [] # will store the position of the errors
    for i in range(0, len(coef_pos)):
        l = 255 - coef_pos[i]
        X.append( gf_pow(2, -l) )

    # Forney algorithm: compute the magnitudes
    E = [0] * (len(msg_in)) # 用来保存误码数据. This is sometimes called the error magnitude polynomial.
    Xlength = len(X)
    for i, Xi in enumerate(X):

        Xi_inv = gf_inverse(Xi)

        # 用errata locator的形式导数(formal derivative)作为Forney算法的分母(denominator), 则第i个error值为： error_evaluator(gf_inverse(Xi)) / error_locator_derivative(gf_inverse(Xi)). 
        err_loc_prime_tmp = []
        for j in range(0, Xlength):
            if j != i:
                err_loc_prime_tmp.append( gf_sub(1, gf_mul(Xi_inv, X[j])) )
        # 下面的乘法计算Forney算法所需的分母 (errata locator derivative)
        err_loc_prime = 1
        for coef in err_loc_prime_tmp:
            err_loc_prime = gf_mul(err_loc_prime, coef)
        # 等价于： err_loc_prime = functools.reduce(gf_mul, err_loc_prime_tmp, 1)

        # Compute y (evaluation of the errata evaluator polynomial)
        # 更直白的Forney算法实现： 
        # Yl = omega(Xl.inverse()) / prod(1 - Xj*Xl.inverse()) for j in len(X)
        y = gf_poly_eval(err_eval[::-1], Xi_inv) # Forney 算法的分子
        y = gf_mul(gf_pow(Xi, 1), y)
        
        # Compute the magnitude
        magnitude = gf_div(y, err_loc_prime) # magnitude value of the error, calculated by the Forney algorithm 
        E[err_pos[i]] = magnitude # store the magnitude for this error into the magnitude polynomial

    # 纠正误码! (纠错码也会被修正!)
    # (这里并不是Forney 算法, 只是直接套用了解码结果)
    msg_in = gf_poly_add(msg_in, E) # 等价于Ci = Ri - Ei ； Ci 为正确信息, Ri 为包含错误的信息, Ei 为修正量表(errata magnitudes) .
    return msg_in

def rs_correct_msg(msg_in, nsym, erase_pos=None):
    '''Reed-Solomon main decoding function'''
    if len(msg_in) > 255: # can't decode, message is too big
        raise ValueError("Message is too long (%i when max is 255)" % len(msg_in))

    msg_out = list(msg_in)     # copy of message
    # erasures: 为了解码是debug方便, 设置为 null (不必要)
    if erase_pos is None:
        erase_pos = []
    else:
        for e_pos in erase_pos:
            msg_out[e_pos] = 0
    # check if there are too many erasures to correct (beyond the Singleton bound)
    if len(erase_pos) > nsym: print("Too many erasures to correct")#raise ReedSolomonError("Too many erasures to correct")
    # prepare the syndrome polynomial using only errors.
    synd = rs_calc_syndromes(msg_out, nsym)
    # 查验 codeword中是否有误码. 如果没有 (伴随式全部系数为 0), codeword原样返回.
    if max(synd) == 0:
        return msg_out[:-nsym], msg_out[-nsym:]  # no errors

    # compute the Forney syndromes, which hide the erasures from the original syndrome (so that BM will just have to deal with errors, not erasures)
    fsynd = rs_forney_syndromes(synd, erase_pos, len(msg_out))
    # compute the error locator polynomial using Berlekamp-Massey
    err_loc = rs_find_error_locator(fsynd, nsym, erase_count=len(erase_pos))
    # locate the message errors using Chien search (or brute-force search)
    err_pos = rs_find_errors(err_loc[::-1] , len(msg_out))
    if err_pos is None:
        #raise ReedSolomonError("Could not locate error")    # error location failed
        print("Could not locate error")    # error location failed

    # Find errors values and apply them to correct the message
    # compute errata evaluator and errata magnitude polynomials, then correct errors and erasures
    msg_out = rs_correct_errata(msg_out, synd, (erase_pos + err_pos)) # note that we here use the original syndrome, not the forney syndrome
                                                                                                                                  # (because we will correct both errors and erasures, so we need the full syndrome)
    # check if the final message is fully repaired
    synd = rs_calc_syndromes(msg_out, nsym)
    if max(synd) > 0:
        #raise ReedSolomonError("Could not correct message")     # message could not be repaired
        print("Could not correct message")     # message could not be repaired
    # return the successfully decoded message
    return msg_out[:-nsym], msg_out[-nsym:] # also return the corrected ecc block so that the user can check()


#############################TESTING#########################################
#############################TESTING#########################################
#############################TESTING#########################################
#msg = [0]*23
#msg_ori = [0,1,0,1,0,0, 1,1, 0,1,0,1,0,0, 1,1,0] #17
#msg_ori = list(range(1,11))+[0,0,0]+list(range(11,15))
msg_ori = list(range(1,18))
init_tables(0x11d)
synd_num = 6 
msg_enc =rs_encode_msg(msg_ori,synd_num)
print("msg_enc is",msg_enc)
#msg_enc_poision = msg_enc[0:1] + [88] + msg_enc[2:10] + [77] + msg_enc[11:19] + [66] + msg_enc[20:23]
#msg_enc_poision = msg_enc[0:1] + [0] + msg_enc[2:10] + [0] + msg_enc[11:19] + [0] + msg_enc[20:23]
#msg_enc_poision = msg_enc[0:1] + [88] + msg_enc[2:10] + [77] + msg_enc[11:23]
msg_enc_poision = msg_enc[0:2] + [88] + msg_enc[3:23]
print("msg_enc_poision is",msg_enc_poision)
syndromes = rs_calc_syndromes(msg_enc_poision,synd_num)
print("syndrome",syndromes)
if max(syndromes) == 0:
    print("ALL RIGHT,DONE\n")
else:
    print("Error Exist\n")

    err_poly = rs_find_error_locator(syndromes,synd_num)
    print("Error PoLy",err_poly[::-1])
    err_pos = rs_find_errors(err_poly[::-1],len(msg_enc))
    print("Error Pos",err_pos)
    print("correct out",rs_correct_errata(msg_enc_poision,syndromes,err_pos))
    #print("correct msg",rs_correct_msg(msg_enc_poision,6))

