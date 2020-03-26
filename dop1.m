clear all;
close all;

crc = [1 0 1 0 0 1 1 0 1 0 1 1 1 1 0 0 1]; %crc16 m=8 x16+x13+x12+x11+x10+x8+x6+x5+x2+1
%crc = [1 0 1 0 1 0 1 1 1]; %x8+x7+x6+x4+x2+1 m=4 crc8
%crc = [1 1 1 0 1]; %x4+x2+x+1 m=3
%crc = [1 1 0 1]; %x3+x+1 m=4
m = 8; %длина информационной последовательности
r = length(crc)-1;
index = 1;
n = 250000;%кол-во тестов
snrArr = -30:5:10;
snrTheorArr = -30:10;
Pe_decodeTheor = ones(1,length(snrTheorArr)).*(1/2^r);%¬еро€тность ошибки декодировани€
Pb_theor = qfunc(sqrt(2*(10.^(snrTheorArr./10))));%¬еро€тность ошибки на бит
codewords = zeros(2^m, m+r);
modulateCodewords = zeros(2^m, m+r);
dArr = zeros((2^m), 1);%массив весов кодовых слов

%таблица кодовых слов
for word = 0:bi2de(ones(1, m))
    [~, c] = gfdeconv(de2bi(bitshift(word, r)), crc);%c(x) = m(x)*(x^r)mod g(x)
    mxr = bitshift(word,r);
    a_dec = bitxor(mxr,bi2de(c));
    a = de2bi(a_dec, m + r);
    codewords(word+1, :) = a;   %кодовое слово
    dArr(word+1) = nnz(a); %добавл€ем слово с ненулевым весом
    a_m = a.*2 - 1;%модул€ци€ кодового слова
    modulateCodewords(word+1, :) = a_m;%добавл€ем смодулированнное кодовое слово в таблицу
end

dMin = min(dArr(2:end));%кодовое слово с наименьшим весом
P_ed = zeros(1, length(Pb_theor));%точна€ веро€тность ошибки декодировани€

%рассчет точной веро€тности ошибки декодировани€
for p = Pb_theor
    tmp = 0;
    for i = dMin:(m+r)% —умма по i = d до n
        Ai = sum(dArr == i);%все кодовые слова с весом i
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i));%промежуточна€ сумма
    end
    P_ed(index) = tmp;
    index = index + 1;
end
index = 1;

%модель канала
for SNR = snrArr
    N_b = 0;%кол-во ошибок на бит
    Ne_decode = 0;%кол-во ошибок декодера
    for i = 1:n
        %создание сообщени€
        rnd_ind = randi(2^m,1);
        massage = modulateCodewords(rnd_ind, :);
        %канал
        SNRi = 10^(SNR/10);%SNR в дЅ
        sigma = sqrt(1/(2*SNRi));
        b_m = massage +sigma*randn(1, length(massage));
        %демодул€тор
        b = b_m > 0;%демодул€ци€ прин€того сообщени€ 
        a = codewords(rnd_ind, :);
        e = xor(a,b);%проверка
        N_b = N_b + nnz(e);%суммируем кол-во ошибок
        [~, s] = gfdeconv(double(b), crc);%синдром
        if (bi2de(s) == 0) & (nnz(e) ~= 0)%ошибка декодера
            Ne_decode = Ne_decode + 1;
        end
    end
    Pb(index) = N_b/(n*(m+r));%¬еро€тность ошибки на бит
    Pe_decode(index) = Ne_decode/n;%¬еро€тность ошибки декодера
    index = index + 1;
end

%веро€тность ошибки на бит
figure(1);
semilogy(snrTheorArr, Pb_theor, 'm', snrArr, Pb, 'ko');
%веро€тноость ошибки декодировани€
figure(2);
semilogy(snrTheorArr, Pe_decodeTheor, 'g', snrArr, Pe_decode, 'ko', snrTheorArr, P_ed, 'm');
