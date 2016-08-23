function dy = BCR_model(t,y,parameter)
dy = zeros(15,1);    

T_syk=1.18;
T_lyn=3.27;


lyn_norm=1.076*exp(-t*0.009744)-0.9535*exp(-t*1.315);
syk_norm=1.309*exp(-t*0.0382)-1.013*exp(-t*0.4757);


pLYN=parameter(1)*T_lyn*lyn_norm;
pSYK=parameter(2)*T_syk*syk_norm;


T_blnk=0.65;
T_bcap=0.9;
T_cd19=0.83;
T_ship=2.82;
T_btk=1.49;
T_pi3k=0.33;
T_pten=0.02;
T_plc=2.57;
T_pip5k=4.4;
T_pkc=2.5;
T_pdk=0.27;
T_akt=0.2;
T_shp1=6.9;
T_ptp1b=1.48;
T_rasgrp=1.03;
T_ras=0.15;
T_raf=1.1;
T_mek=1.39;
T_erk=2.26;
T_d=1.0;
T_pip=5;


k_bcap=parameter(3);K_bcap=parameter(28);
k_shbl=parameter(4);K_shbl=parameter(29);
k_ptsh=parameter(5);K_ptsh=parameter(30);
k_ptpl=parameter(6);K_ptpl=parameter(31);
k_dbtk=parameter(7);K_dbtk=parameter(32);
k_cd19=parameter(8);K_cd19=parameter(33);
k_dakt=parameter(9);K_dakt=parameter(34);
k_dmek1=parameter(10);K_dmek1=parameter(35);
k_dmek2=parameter(11);K_dmek2=parameter(36);
k_drasgrp=parameter(12);K_drasgrp=parameter(37);
k_derk=parameter(13);K_derk=parameter(38);

k_lyn=parameter(14); K_lyn=parameter(39);
k_syk=parameter(15); K_syk=parameter(40);
k_btk=parameter(16); K_btk=parameter(41);
k_pipi=parameter(17); K_pipi=parameter(42);
k_ptpi=parameter(18); K_ptpi=parameter(43);
k_pi3k=parameter(19); K_pi3k=parameter(44);
k_plpi=parameter(20); K_plpi=parameter(45);
k_shpi=parameter(21); K_shpi=parameter(46);
k_pdak=parameter(22); K_pdak=parameter(47);
k_pkc=parameter(23); K_pkc=parameter(48);
k_rara=parameter(24);K_rara=parameter(49);
k_rame=parameter(25);K_rame=parameter(50);
k_erme=parameter(26);K_erme=parameter(51);
k_meer=parameter(27);K_meer=parameter(52);

k_b4=parameter(53);
k_b5=parameter(54);
k_b6=parameter(55);
k_mept=parameter(56);
k_b1=parameter(57); k_b2=parameter(58); 
k_b3=parameter(59); 

k_dapk=parameter(60);
k_pdk=parameter(61);
k_akt=parameter(62);
k_dara=parameter(63);
k_dras=parameter(64);
k_raf=parameter(65);


g_dag=parameter(66)
k=parameter(67);
g_pip2=parameter(68);
K_erra=parameter(69);K_akra=parameter(70);K_pkbt=parameter(71);%1;



pBLNK=y(1);
pBCAP=y(2);
pCD19=y(3);
pSHIP=y(4);
pBTK=y(5);
pPLC=y(6);
PIP2=y(7);
PIP3=y(8);
DAG=y(9);
pAKT=y(10);
pRASGRP=y(11);
RAS=y(12);
pMEK=y(13);
pMEK_=y(14);
pERK=y(15);


PTEN=(1+k_mept*(pMEK_+T_pten)-sqrt((1+k_mept*(pMEK_+T_pten))^2-4*k_mept^2*pMEK_*T_pten))/(2*k_mept);  
PKC=(1+k_dapk*(DAG+T_pkc)-sqrt((1+k_dapk*(DAG+T_pkc))^2-4*k_dapk^2*DAG*T_pkc))/(2*k_dapk); 
BTK=((1+k_b1*(pBLNK+T_btk)-sqrt((1+k_b1*(pBLNK+T_btk))^2-4*k_b1^2*pBLNK*T_btk))/(2*k_b1)+(1+k_b2*(PIP3+T_btk)-sqrt((1+k_b2*(PIP3+T_btk))^2-4*k_b2^2*PIP3*T_btk))/(2*k_b2))/(1+K_pkbt*PKC);
PLC=(1+k_b6*(pBLNK+T_plc)-sqrt((1+k_b6*(pBLNK+T_plc))^2-4*k_b6^2*pBLNK*T_plc))/(2*k_b6);
PI3K=(1+k_b3*(pBCAP+T_pi3k)-sqrt((1+k_b3*(pBCAP+T_pi3k))^2-4*k_b3^2*pBCAP*T_pi3k))/(2*k_b3)+(1+k_b4*(pCD19+T_pi3k)-sqrt((1+k_b4*(pCD19+T_pi3k))^2-4*k_b4^2*pCD19*T_pi3k))/(2*k_b4);
PIP5K=(1+k_b5*(BTK+T_pip5k)-sqrt((1+k_b5*(BTK+T_pip5k))^2-4*k_b5^2*BTK*T_pip5k))/(2*k_b5);
PDK=(1+k_pdk*(PIP3+T_pdk)-sqrt((1+k_pdk*(PIP3+T_pdk))^2-4*k_pdk^2*PIP3*T_pdk))/(2*k_pdk);
AKT=(1+k_akt*(PIP3+T_akt)-sqrt((1+k_akt*(PIP3+T_akt))^2-4*k_akt^2*PIP3*T_akt))/(2*k_akt);
RASGRP=(1+k_dara*(DAG+T_rasgrp)-sqrt((1+k_dara*(DAG+T_rasgrp))^2-4*k_dara^2*DAG*T_rasgrp))/(2*k_dara);
RAF=(1+k_raf*(RAS+T_raf)-sqrt((1+k_raf*(RAS+T_raf))^2-4*k_raf^2*RAS*T_raf))/(2*k_raf)/(1+K_erra*pERK+K_akra*pAKT);


dy(1) = k_syk*pSYK*(T_blnk-pBLNK)/(K_syk+T_blnk-pBLNK)-k_shbl*T_shp1*pBLNK/(K_shbl+pBLNK);
dy(2) = k_syk*pSYK*(T_bcap-pBCAP)/(K_syk+T_bcap-pBCAP)+k_btk*pBTK*(T_bcap-pBCAP)/(K_btk+T_bcap-pBCAP)-k_bcap*T_d*pBCAP/(K_bcap+pBCAP);
dy(3) = k_lyn*pLYN*(T_cd19-pCD19)/(K_lyn+T_cd19-pCD19)-k_cd19*T_d*pCD19/(K_cd19+pCD19);
dy(4) = k_lyn*pLYN*(T_ship-pSHIP)/(K_lyn+T_ship-pSHIP)-k_ptsh*T_ptp1b*pSHIP/(K_ptsh+pSHIP);
dy(5) = k_syk*pSYK*(BTK-pBTK)/(K_syk+BTK-pBTK)+k_lyn*pLYN*(BTK-pBTK)/(K_lyn+BTK-pBTK)-k_dbtk*T_d*pBTK/(K_dbtk+pBTK);
dy(6) = k_btk*pBTK*(PLC-pPLC)/(K_btk+PLC-pPLC)-k_ptpl*T_ptp1b*pPLC/(K_ptpl+pPLC);
dy(7) = k-g_pip2*PIP2+k_pipi*PIP5K*T_pip/(K_pipi+T_pip)+k_ptpi*PTEN*PIP3/(K_ptpi+PIP3)-k_pi3k*PI3K*PIP2/(K_pi3k+PIP2)-k_plpi*pPLC*PIP2/(K_plpi+PIP2);
dy(8) = k_pi3k*PI3K*PIP2/(K_pi3k+PIP2)-k_ptpi*PTEN*PIP3/(K_ptpi+PIP3)-k_shpi*pSHIP*PIP3/(K_shpi+PIP3);
dy(9) = k_plpi*pPLC*PIP2/(K_plpi+PIP2)-g_dag*DAG;
dy(10)= k_pdak*PDK*(AKT-pAKT)/(K_pdak+AKT-pAKT)-k_dakt*T_d*pAKT/(K_dakt+pAKT);
dy(11)= k_pkc*PKC*(RASGRP-pRASGRP)/(K_pkc+RASGRP-pRASGRP)-k_drasgrp*T_d*pRASGRP/(K_drasgrp+pRASGRP);
dy(12)= k_rara*pRASGRP*(T_ras-RAS)/(K_rara+T_ras-RAS)-k_dras*RAS;
dy(13)= k_rame*RAF*(T_mek-pMEK)/(K_rame+T_mek-pMEK)-k_dmek1*T_d*pMEK/(K_dmek1+pMEK);
dy(14)= k_erme*pERK*(T_mek-pMEK_)/(K_erme+T_mek-pMEK_)-k_dmek2*T_d*pMEK_/(K_dmek2+pMEK_);
dy(15)= k_meer*pMEK*(T_erk-pERK)/(K_meer+T_erk-pERK)-k_derk*T_d*pERK/(K_derk+pERK);



