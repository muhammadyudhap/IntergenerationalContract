function f = dW0(XX,parameters);

parameters0(1)      = parameters(1);
parameters0(2)      = parameters(2);
parameters0(3)      = parameters(3);
parameters0(4)      = parameters(4);
parameters0(5)      = parameters(5);
parameters0(6)      = 0;
parameters0(7)      = parameters(7);
parameters0(8)      = parameters(8);
parameters0(9)      = parameters(9);
parameters0(10)     = parameters(10);
parameters0(11)     = parameters(11);
parameters0(12)     = parameters(12);
parameters0(13)     = parameters(13);
parameters0(14)     = parameters(14);
parameters0(15)     = parameters(15);

parameters1(1)      = parameters(1);
parameters1(2)      = parameters(2);
parameters1(3)      = parameters(3);
parameters1(4)      = parameters(4);
parameters1(5)      = parameters(5);
parameters1(6)      = parameters(6);
parameters1(7)      = parameters(7);
parameters1(8)      = parameters(8);
parameters1(9)      = XX;
parameters1(10)     = parameters(10);
parameters1(11)     = parameters(11);
parameters1(12)     = parameters(12);
parameters1(13)     = parameters(13);
parameters1(14)     = parameters(14);
parameters1(15)     = parameters(15);

[U0A, U1A] = welfarelog(parameters0);
[U0B, U1B] = welfarelog(parameters1);

f = U1A-U1B;

