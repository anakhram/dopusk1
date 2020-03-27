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
Pe_decodeTheor = ones(1,length(snrTheorArr)).*(1/2^r);%Вероятность ошибки декодирования
Pb_theor = qfunc(sqrt(2*(10.^(snrTheorArr./10))));%Вероятность ошибки на бит
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
    dArr(word+1) = nnz(a); %добавляем слово с ненулевым весом
    a_m = a.*2 - 1;%модуляция кодового слова
    modulateCodewords(word+1, :) = a_m;%добавляем смодулированнное кодовое слово в таблицу
end

dMin = min(dArr(2:end));%кодовое слово с наименьшим весом
P_ed = zeros(1, length(Pb_theor));%точная вероятность ошибки декодирования

%рассчет точной вероятности ошибки декодирования
for p = Pb_theor
    tmp = 0;
    for i = dMin:(m+r)% Сумма по i = d до n
        Ai = sum(dArr == i);%все кодовые слова с весом i
        tmp = tmp + (Ai * p^i * (1-p)^((m+r)- i));%промежуточная сумма
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
        %создание сообщения
        rnd_ind = randi(2^m,1);
        massage = modulateCodewords(rnd_ind, :);
        %канал
        SNRi = 10^(SNR/10);%SNR в дБ
        sigma = sqrt(1/(2*SNRi));
        b_m = massage +sigma*randn(1, length(massage));
        %демодулятор
        b = b_m > 0;%демодуляция принятого сообщения 
        e = xor(codewords(rnd_ind, :),b);%проверка
        N_b = N_b + nnz(e);%суммируем кол-во ошибок
        [~, s] = gfdeconv(double(b), crc);%синдром
        if (bi2de(s) == 0) & (nnz(e) ~= 0)%ошибка декодера
            Ne_decode = Ne_decode + 1;
        end
    end
    Pb(index) = N_b/(n*(m+r));%Вероятность ошибки на бит
    Pe_decode(index) = Ne_decode/n;%Вероятность ошибки декодера
    index = index + 1;
end

%вероятность ошибки на бит
figure(1);
semilogy(snrTheorArr, Pb_theor, 'm', snrArr, Pb, 'ko');
%вероятноость ошибки декодирования
figure(2);
semilogy(snrTheorArr, Pe_decodeTheor, 'g', snrArr, Pe_decode, 'ko', snrTheorArr, P_ed, 'm');
