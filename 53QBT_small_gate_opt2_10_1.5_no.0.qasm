OPENQASM 2.0;
include "qelib1.inc";
qreg q[53];
cx q[26],q[48];
cx q[44],q[24];
h q[12];
h q[5];
cx q[5],q[12];
cx q[31],q[0];
h q[24];
h q[12];
cx q[33],q[20];
h q[5];
h q[42];
h q[26];
cx q[17],q[42];
h q[28];
cx q[31],q[0];
h q[2];
h q[26];
h q[48];
h q[31];
h q[16];
cx q[2],q[16];
h q[20];
cx q[2],q[28];
cx q[16],q[28];
h q[20];
cx q[20],q[16];
h q[48];
h q[20];
cx q[41],q[36];
h q[48];
h q[17];
cx q[15],q[31];
cx q[10],q[50];
cx q[26],q[48];
h q[41];
cx q[38],q[8];
cx q[20],q[16];
h q[28];
h q[15];
h q[16];
h q[10];
h q[26];
cx q[17],q[42];
h q[50];
h q[10];
h q[50];
h q[44];
h q[33];
h q[28];
h q[25];
cx q[20],q[2];
cx q[24],q[31];
h q[34];
cx q[20],q[28];
h q[20];
h q[31];
cx q[17],q[42];
h q[21];
cx q[25],q[38];
h q[3];
cx q[50],q[34];
h q[38];
h q[33];
cx q[39],q[44];
h q[34];
h q[25];
h q[2];
cx q[50],q[3];
cx q[21],q[33];
h q[24];
cx q[50],q[3];
h q[52];
cx q[39],q[31];
h q[24];
h q[39];
cx q[7],q[23];
cx q[49],q[42];
h q[7];
h q[52];
cx q[29],q[19];
h q[42];
h q[24];
h q[24];
cx q[24],q[39];
cx q[42],q[52];
h q[23];
h q[7];
h q[52];
h q[39];
cx q[24],q[39];
h q[39];
h q[49];
cx q[24],q[39];
cx q[23],q[43];
cx q[49],q[52];
h q[52];
h q[12];
cx q[24],q[44];
h q[27];
h q[44];
h q[24];
h q[31];
h q[31];
h q[5];
h q[20];
cx q[5],q[12];
h q[31];
cx q[5],q[12];
cx q[27],q[25];
h q[1];
cx q[20],q[2];
h q[2];
cx q[20],q[2];
h q[30];
cx q[30],q[1];
cx q[31],q[15];
h q[31];
h q[25];
cx q[0],q[52];
h q[19];
cx q[36],q[31];
cx q[27],q[25];
h q[36];
h q[3];
h q[46];
h q[30];
h q[0];
cx q[49],q[0];
h q[52];
cx q[36],q[31];
h q[0];
cx q[32],q[30];
h q[0];
cx q[18],q[46];
h q[46];
h q[46];
h q[46];
h q[25];
cx q[1],q[19];
cx q[32],q[30];
h q[31];
cx q[10],q[3];
