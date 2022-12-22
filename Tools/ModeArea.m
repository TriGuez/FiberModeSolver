function Aeff = ModeArea(x,y,Field)

Aeff = (trapz(x,trapz(y,abs(Field).^2)).^2)./((trapz(x,trapz(y,abs(Field).^4))));