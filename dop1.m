clear all;
close all;

crc = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1]; %crc16 m=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
%crc = [1 0 1 0 1 0 1 1 1]; %x8+x7+x6+x4+x2+1 m=4 crc8
%crc = [1 1 1 0 1]; %x4+x2+x+1 m=3
%crc = [1 1 0 1]; %x3+x+1 m=4
m = 8; %����� �������������� ������������������
r = length(crc)-1;
index = 1;
n = 250000;%���-�� ������
snrArr = -30:5:10;
snrTheorArr = -30:10;
Pe_decodeTheor = ones(1,length(snrTheorArr)).*(1/2^r);%����������� ������ �������������
Pb_theor = qfunc(sqrt(2*(10.^(snrTheorArr./10))));%����������� ������ �� ���
codewords = zeros(2^m, m+r);
modulateCodewords = zeros(2^m, m+r);
dArr = zeros((2^m), 1);%������ ����� ������� ����

%������� ������� ����
for word = 0:bi2de(ones(1, m))
    [~, c] = gfdeconv(de2bi(bitshift(word, r)), crc);%c(x) = m(x)*(x^r)mod g(x)
    mxr = bitshift(word,r);
    a_dec = bitxor(mxr,bi2de(c));
    a = de2bi(a_dec, m + r);
    codewords(word+1, :) = a;   %������� �����
    dArr(word+1) = nnz(a); %��������� ����� � ��������� �����
    a_m = a.*2 - 1;%��������� �������� �����
    modulateCodewords(word+1, :) = a_m;%��������� ���������������� ������� ����� � �������
end

dMin = min(dArr(2:end));%������� ����� � ���������� �����
P_ed = zeros(1, length(Pb_theor));%������ ����������� ������ �������������

%������� ������ ����������� ������ �������������
for p = Pb_theor
    tmp = 0;
    for i = dMin:(m+r)% ����� �� i = d �� n
        Ai = sum(dArr == i);%��� ������� ����� � ����� i
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i));%������������� �����
    end
    P_ed(index) = tmp;
    index = index + 1;
end
index = 1;

%������ ������
for SNR = snrArr
    N_b = 0;%���-�� ������ �� ���
    Ne_decode = 0;%���-�� ������ ��������
    for i = 1:n
        %�������� ���������
        rnd_ind = randi(2^m,1);
        massage = modulateCodewords(rnd_ind, :);
        %�����
        SNRi = 10^(SNR/10);%SNR � ��
        sigma = sqrt(1/(2*SNRi));
        b_m = massage +sigma*randn(1, length(massage));
        %�����������
        b = b_m > 0;%����������� ��������� ��������� 
        a = codewords(rnd_ind, :);
        e = xor(a,b);%��������
        N_b = N_b + nnz(e);%��������� ���-�� ������
        [~, s] = gfdeconv(double(b), crc);%�������
        if (bi2de(s) == 0) & (nnz(e) ~= 0)%������ ��������
            Ne_decode = Ne_decode + 1;
        end
    end
    Pb(index) = N_b/(n*(m+r));%����������� ������ �� ���
    Pe_decode(index) = Ne_decode/n;%����������� ������ ��������
    index = index + 1;
end

%����������� ������ �� ���
figure(1);
semilogy(snrTheorArr, Pb_theor, 'm', snrArr, Pb, 'ko');
%������������ ������ �������������
figure(2);
semilogy(snrTheorArr, Pe_decodeTheor, 'g', snrArr, Pe_decode, 'ko', snrTheorArr, P_ed, 'm');
