function [f] = inverseConvOld(f,g)
y = conv(f,g);
n = length(y);
nn = length(f);
r = xcorr(f,g);

ff = inverseConv(g,y);

i = 1;

figure(5)

subplot(4,4,i)
i = i+1;
plot(y)
title("conv(f,g)")

subplot(4,4,i)
i = i+1;
plot(r)
title("xcorr(f,g)")

subplot(4,4,i)
i = i+1;
plot(f)
title("f and ff")
hold on
plot(ff(1:nn))
hold off

subplot(4,4,i)
i = i+1;
plot(g)
title("g")

subplot(4,4,i)
i = i+1;
plot(ff)
title("ff")

