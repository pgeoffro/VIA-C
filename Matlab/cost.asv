function [C,gradC,HessC] = cost()

syms a b c Winst xdes ydes

C = Winst*( (cos(a) + cos(a+b) + cos(a+b+c) - xdes)^2 + (sin(a) + sin(a+b) + sin(a+b+c) - ydes)^2);

gradC = [diff(C,a);diff(C,b);diff(C,c)];

HessC = [diff(diff(C,a),a),diff(diff(C,a),b),diff(diff(C,a),c);diff(diff(C,b),a),diff(diff(C,b),b),diff(diff(C,c),b);diff(diff(C,a),c),diff(diff(C,b),c),diff(diff(C,c),c)];