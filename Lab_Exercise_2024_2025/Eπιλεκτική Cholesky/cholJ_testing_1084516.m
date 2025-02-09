n = 100;
rng('default'); % Ή rng(id) όπου id είναι το seed
d = rand(n,1);  % Διανύσμα d με θετικές τιμές
u = randn(n,1);

times = zeros(n,1);
for j = 1:n
    tic;
    [ljj, ljj_1] = cholJ_1084516(u,d,j);
    times(j) = toc;
end

figure;
plot(1:n, times, 'LineWidth',2);
xlabel('j');
ylabel('Χρόνος εκτέλεσης (sec)');
title('Χρονισμός της cholJ\_ID για j = 1:n');
