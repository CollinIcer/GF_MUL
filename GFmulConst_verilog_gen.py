from bitarray import bitarray
P=[1,0,1,1,1,0,0,0,1]     #primitive : X^8+X^4+X^3+X^2+1
Snum = 8 #simbol bitwidth 
tmp=[0]*Snum

def GFmulConst_verilog_gen(Const,P,fp):
    Cbin = bin(Const)
    Cbin = bitarray(Cbin[2:])
    Pbin = bitarray(P[:-1])
    bitarray.reverse(Cbin) 

    while(len(Cbin)<Snum):
        Cbin.append(0)

    for i in range(Snum): 
        if(i==0):
            tmp[i] = Cbin
        else:
            tmp[i] = tmp[i-1][0:Snum-1]
            tmp[i].insert(0,0)
            if(tmp[i-1][7]):
                tmp[i] = tmp[i] ^ Pbin

    #print(tmp)

    fp.write('module gf_mul%d( input [%d:0] din, output [%d:0] dout);\n'%(Const,Snum-1,Snum-1))
    for i in range(Snum):
        fp.write('assign dout[%d] = 0 '%i)
        for j in range(Snum):
            if(tmp[j][i]):
                fp.write(' ^ din[%d]'%j)
        fp.write(';\n')
    fp.write('endmodule\n')
        
fp = open("gf_mul.v",'w')
GFmulConst_verilog_gen(126,P,fp) 
GFmulConst_verilog_gen(4,P,fp) 
GFmulConst_verilog_gen(158,P,fp) 
GFmulConst_verilog_gen(28,P,fp) 
GFmulConst_verilog_gen(49,P,fp) 
GFmulConst_verilog_gen(117,P,fp) 
fp.close()