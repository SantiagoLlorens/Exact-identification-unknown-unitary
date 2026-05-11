function o=odskok(b,a)
global D d n 
bitov=round(log2(b))+1;
suma=0;
base=2^(2*a);
for k = 1:bitov
    suma=suma+(2^(2*(k-1)))*bitget(b-1,k);
end;
o=base*suma;