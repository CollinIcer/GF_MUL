module gf_mul126( input [7:0] din, output [7:0] dout);
assign dout[0] = 0  ^ din[2] ^ din[3] ^ din[4] ^ din[5] ^ din[7];
assign dout[1] = 0  ^ din[0] ^ din[3] ^ din[4] ^ din[5] ^ din[6];
assign dout[2] = 0  ^ din[0] ^ din[1] ^ din[2] ^ din[3] ^ din[6];
assign dout[3] = 0  ^ din[0] ^ din[1] ^ din[5];
assign dout[4] = 0  ^ din[0] ^ din[1] ^ din[3] ^ din[4] ^ din[5] ^ din[6] ^ din[7];
assign dout[5] = 0  ^ din[0] ^ din[1] ^ din[2] ^ din[4] ^ din[5] ^ din[6] ^ din[7];
assign dout[6] = 0  ^ din[0] ^ din[1] ^ din[2] ^ din[3] ^ din[5] ^ din[6] ^ din[7];
assign dout[7] = 0  ^ din[1] ^ din[2] ^ din[3] ^ din[4] ^ din[6] ^ din[7];
endmodule
module gf_mul4( input [7:0] din, output [7:0] dout);
assign dout[0] = 0  ^ din[6];
assign dout[1] = 0  ^ din[7];
assign dout[2] = 0  ^ din[0] ^ din[6];
assign dout[3] = 0  ^ din[1] ^ din[6] ^ din[7];
assign dout[4] = 0  ^ din[2] ^ din[6] ^ din[7];
assign dout[5] = 0  ^ din[3] ^ din[7];
assign dout[6] = 0  ^ din[4];
assign dout[7] = 0  ^ din[5];
endmodule
module gf_mul158( input [7:0] din, output [7:0] dout);
assign dout[0] = 0  ^ din[1] ^ din[4];
assign dout[1] = 0  ^ din[0] ^ din[2] ^ din[5];
assign dout[2] = 0  ^ din[0] ^ din[3] ^ din[4] ^ din[6];
assign dout[3] = 0  ^ din[0] ^ din[5] ^ din[7];
assign dout[4] = 0  ^ din[0] ^ din[4] ^ din[6];
assign dout[5] = 0  ^ din[1] ^ din[5] ^ din[7];
assign dout[6] = 0  ^ din[2] ^ din[6];
assign dout[7] = 0  ^ din[0] ^ din[3] ^ din[7];
endmodule
module gf_mul28( input [7:0] din, output [7:0] dout);
assign dout[0] = 0  ^ din[4] ^ din[5] ^ din[6];
assign dout[1] = 0  ^ din[5] ^ din[6] ^ din[7];
assign dout[2] = 0  ^ din[0] ^ din[4] ^ din[5] ^ din[7];
assign dout[3] = 0  ^ din[0] ^ din[1] ^ din[4];
assign dout[4] = 0  ^ din[0] ^ din[1] ^ din[2] ^ din[4] ^ din[6];
assign dout[5] = 0  ^ din[1] ^ din[2] ^ din[3] ^ din[5] ^ din[7];
assign dout[6] = 0  ^ din[2] ^ din[3] ^ din[4] ^ din[6];
assign dout[7] = 0  ^ din[3] ^ din[4] ^ din[5] ^ din[7];
endmodule
module gf_mul49( input [7:0] din, output [7:0] dout);
assign dout[0] = 0  ^ din[0] ^ din[3] ^ din[4] ^ din[7];
assign dout[1] = 0  ^ din[1] ^ din[4] ^ din[5];
assign dout[2] = 0  ^ din[2] ^ din[3] ^ din[4] ^ din[5] ^ din[6] ^ din[7];
assign dout[3] = 0  ^ din[5] ^ din[6];
assign dout[4] = 0  ^ din[0] ^ din[3] ^ din[4] ^ din[6];
assign dout[5] = 0  ^ din[0] ^ din[1] ^ din[4] ^ din[5] ^ din[7];
assign dout[6] = 0  ^ din[1] ^ din[2] ^ din[5] ^ din[6];
assign dout[7] = 0  ^ din[2] ^ din[3] ^ din[6] ^ din[7];
endmodule
module gf_mul117( input [7:0] din, output [7:0] dout);
assign dout[0] = 0  ^ din[0] ^ din[2] ^ din[3] ^ din[4];
assign dout[1] = 0  ^ din[1] ^ din[3] ^ din[4] ^ din[5];
assign dout[2] = 0  ^ din[0] ^ din[3] ^ din[5] ^ din[6];
assign dout[3] = 0  ^ din[1] ^ din[2] ^ din[3] ^ din[6] ^ din[7];
assign dout[4] = 0  ^ din[0] ^ din[7];
assign dout[5] = 0  ^ din[0] ^ din[1];
assign dout[6] = 0  ^ din[0] ^ din[1] ^ din[2];
assign dout[7] = 0  ^ din[1] ^ din[2] ^ din[3];
endmodule