function h_output = h(theta)
   gammaf = theta(1);
   lambda = theta(2);
   rhopi = theta(3);
   rhox = theta(4);

   h1 = (rhox - rhox*rhopi*gammaf)/(rhox*gammaf + lambda);
   h2 = rhox;
   h3 = (rhopi - rhopi^2*gammaf + gammaf - 1)/(rhox*gammaf + lambda);
   h4 = rhopi;

   h_output = [h1 h2 h3 h4]';
end