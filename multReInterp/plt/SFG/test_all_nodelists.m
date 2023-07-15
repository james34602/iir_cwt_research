clc
H01=flow_tf('biquad.flw');
H01.sym{1}
fprintf('\n\nPress space to continue')
pause,clc
H02=flow_tf('CRFBSYM.flw');
H02.sym{1}
H02.sym{2}
fprintf('\n\nPress space to continue')
pause,clc

H03=flow_tf('iir.flw');
H03.sym{1}
fprintf('\n\nPress space to continue')
pause,clc

H04=flow_tf('MODL.flw');
H04.sym{1}
H04.sym{2}
fprintf('\n\nPress space to continue')
pause,clc

H05=flow_tf('CRFBNUM.flw');
H05.sym{1}
H05.sym{2}
H05.tf{1}
H05.tf{2}
fprintf('\n\nPress space to continue')
pause,clc

H06=flow_tf('efb.flw');
H06.sym{1}
H06.sym{2}
H06.tf{1}
H06.tf{2}
fprintf('\n\nPress space to continue')
pause,clc

H07=flow_tf('integra.flw');
H07.sym{1}
H07.tf{1}
fprintf('\n\nPress space to continue')
pause,clc

H08=flow_tf('fiori.flw');
H08.sym{1}
H08.sym{2}
H08.tf{1}
H08.tf{2}
fprintf('\n\nPress space to continue')
pause,clc

H09=flow_tf('MOD1.flw');
H09.sym{1}
H09.sym{2}
H09.tf{1}
H09.tf{2}
fprintf('\n\nPress space to continue')
pause,clc

H10=flow_tf('zpar.flw');
H10.sym{1}
H10.tf{1}
