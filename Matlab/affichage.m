
adresseRT='C:\Users\Perle\Desktop\VIA\RTorque.txt';
RTorque = textread(adresseRT);
plot(RTorque,'g');
hold on;

adresseHT='C:\Users\Perle\Desktop\VIA\HTorque.txt';
HTorque = textread(adresseHT);
plot(HTorque,'k');