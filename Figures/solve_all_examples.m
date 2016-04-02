%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Repeated moral hazard and recursive Lagrangeans
%
%   This script solves all the examples in section 5
%   and generates all figures included in the paper
%
%                       Antonio Mele
%                       March 2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% standard repeated moral hazard model
solve_RMH;
print(1,'-depsc2', 'figure1.eps');
print(2,'-depsc2', 'figure2.eps');
print(3,'-depsc2', 'figure3.eps');
print(4,'-depsc2', 'figure4.eps');
print(5,'-depsc2', 'figure5.eps');



% repeated moral hazard with hidden assets
solve_HA;
print(6,'-depsc2', 'figure6.eps');
print(7,'-depsc2', 'figure7.eps');
print(8,'-depsc2', 'figure8.eps');
print(9,'-depsc2', 'figure9.eps');
print(10,'-depsc2', 'figure10.eps');


% dynamic risk sharing with two-sided moral
% hazard in an endowment economy
solve_RSE;

print(11,'-depsc2', 'figure11.eps');
print(12,'-depsc2', 'figure12.eps');
print(13,'-depsc2', 'figure13.eps');
print(14,'-depsc2', 'figure14.eps');
print(15,'-depsc2', 'figure15.eps');



% dynamic risk sharing with two-sided moral hazard
%in a production economy with capital
solve_RSP;
print(16,'-depsc2', 'figure16.eps');
print(17,'-depsc2', 'figure17.eps');
print(18,'-depsc2', 'figure18.eps');
print(19,'-depsc2', 'figure19.eps');


