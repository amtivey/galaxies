   10 REM Galaxy Collision
   20 REM by Andrew Rankin
   30 REM for B+/M/C/A
   40 REM (c) BAU July 1989
   50 :
   60 MODE 128
   70 v$="1.02"
   80 DIM x(750),y(750),z(750),u(750),v(750),w(750)
   90 galr=30:galm=5:galri=10:galc%=12
  100 mz%=2
  110 phase%=1:@%=&50A
  120 title$="COLLIDING GALAXIES SIMULATION version "+v$
  130 REPEAT
  140 ON phase% PROCsetup,PROCsimulation ELSE
  150 UNTIL phase%=3
  160 END
  170 :
  180 DEF PROCsetup
  190 CLS
  200 VDU 19,1,4;0;19,0,7;0;
  210 PRINT TAB(19,0)title$
  220 PRINT '" This program simulates the passing of two galaxies in close proximity. To"
  230 PRINT " simplify calculations of star positions it is assumed that all of the"
  240 PRINT " galactic mass is at the centre and ignores the gravitational field due to the"
  250 PRINT " surrounding stars. It is assumed the effects of forces other than gravity"
  260 PRINT " are negligible. The target galaxy is taken to be a disk of stars in stable"
  270 PRINT " circular orbit with zero drift velocity and the intruder is taken to be a"
  280 PRINT " point of comparable mass."
  290 :
  300 PRINT TAB(1,12)"SET UP:"TAB(1,13)"In the target galaxy:"
  310 INPUT TAB(10,14)"How many rings of stars are required ",rings%
  320 INPUT TAB(10,15)"How many stars per ring are required ",stpr%
  330 stars%=stpr%*rings%:dr=(galr-galri)/(rings%-1)
  340 PRINT " The target galaxy is located at coordinates (0,0,0) with zero";
  350 PRINT " drift velocity,"
  360 PRINT " and mass 5.0 mass units. ";stars%" stars rotate in circular ";
  370 PRINT "orbits of 10 to 30 distance units."
  380 PRINT TAB(1,19)"For the intruder galaxy:"
  390 INPUT TAB(10,20)"Mass as a fraction (%) of target galaxy ",massfr
  400 INPUT TAB(10,21)"Distance from target galaxy:"TAB(10,22)"X (-70 to 70) ? "d
  410 INPUT TAB(33,22)"Y (-70 to 70) ? "e TAB(56,22)"Z (-35 to 35) ? "f
  420 INPUT TAB(10,23)"Velocity components:"TAB(10,24)"X ? "r TAB(33,24)"Y ? "s
  430 INPUT TAB(56,24)"Z ? "t
  440 PRINT TAB(5,26)"Set up complete - select (A)lter set up, (R)un simulation";
  450 PRINT ", (E)xit [ ]"TAB(74,26);
  460 REPEAT
  470 ans$=GET$
  480 PRINT TAB(74,26)ans$;
  490 phase%=INSTR("ARE",ans$)
  500 UNTIL phase%>0
  510 stop%=-1
  520 ENDPROC
  530 :
  540 DEF PROCsimulation
  550 VDU 12,20,19,0,4;0;23,1,0,0,0,0,0,0,0,0
  560 m=galm
  570 n=m*massfr/100
  580 a=150:b=100:c=0:d=a+d:e=e+b:f=f+c
  590 o=0:p=0:q=0
  600 T%=0:I%=0:sf=2
  610 PRINT TAB(19,0)title$TAB(15,1)"     Target galaxy     "
  620 PRINT TAB(49,1)"     Intruder galaxy     "
  630 PRINT TAB(17,2)"- start:"TAB(35,2)"- now:"
  640 PRINT TAB(51,2)"- start:"TAB(69,2)"- now:"
  650 PRINT TAB(1,3)"Mass"TAB(15,3)m TAB(31,3)m TAB(49,3)n TAB(65,3)n
  660 PRINT TAB(1,4)"X coordinate"TAB(1,5)"Y coordinate"TAB(1,6)"Z coordinate"
  670 PRINT TAB(1,7)"Velocity (X)"TAB(1,8)"Velocity (Y)"TAB(1,9)"Velocity (Z)"
  680 PRINT TAB(1,11)"X-Y plane"TAB(41,11)"X-Z plane (Z shown x 2)"
  690 PRINT TAB(70,11)"Time"TAB(15,4)a TAB(49,4)d TAB(15,5)b TAB(49,5)e
  700 PRINT TAB(15,6)c TAB(49,6)f TAB(15,7)o TAB(49,7)r
  710 PRINT TAB(15,8)p TAB(49,8)s TAB(15,9)q TAB(49,9)t TAB(75,11);T%
  720 FOR J%=0 TO rings%-1
  730 rad=J%*dr+galri
  740 vel=SQR(m/rad)
  750 theta=0.5*vel/rad
  760 IF J%=0 THEN vel=0.9*vel
  770 FOR K%=0 TO stpr%-1
  780 g=K%*2*PI/stpr%
  790 x(I%)=rad*COS(g)+150
  800 y(I%)=rad*SIN(g)+100
  810 z(I%)=0
  820 w(I%)=0
  830 v(I%)=vel*COS(g-theta)
  840 u(I%)=-vel*SIN(g-theta)
  850 I%=I%+1
  860 NEXT:NEXT
  870 PROCplotgalaxies
  880 PROCkey
  890 REPEAT
  900 FOR I%=0 TO stars%-1
  910 f1=m/((x(I%)-a)^2+(y(I%)-b)^2+(z(I%)-c)^2+sf)^1.5
  920 f2=n/((x(I%)-d)^2+(y(I%)-e)^2+(z(I%)-f)^2+sf)^1.5
  930 ax=f1*(a-x(I%))+f2*(d-x(I%))
  940 ay=f1*(b-y(I%))+f2*(e-y(I%))
  950 az=f1*(c-z(I%))+f2*(f-z(I%))
  960 u(I%)=u(I%)+ax
  970 v(I%)=v(I%)+ay
  980 w(I%)=w(I%)+az
  990 x(I%)=x(I%)+u(I%)
 1000 y(I%)=y(I%)+v(I%)
 1010 z(I%)=z(I%)+w(I%)
 1020 NEXT
 1030 rad=((a-d)^2+(b-e)^2+(c-f)^2+sf)^1.5
 1040 ax=(d-a)/rad
 1050 ay=(e-b)/rad
 1060 az=(f-c)/rad
 1070 o=n*ax+o
 1080 p=n*ay+p
 1090 q=n*az+q
 1100 r=r-m*ax
 1110 s=s-m*ay
 1120 t=t-m*az
 1130 a=a+o
 1140 b=b+p
 1150 c=c+q
 1160 d=d+r
 1170 e=e+s
 1180 f=f+t
 1190 T%=T%+1
 1200 PROCplotgalaxies
 1210 key%=INKEY(5)
 1220 IF key%=80 OR T%=stop% PROCwait
 1230 IF key%=73 PROCiterate
 1240 IF key%=71 PROCsdump:key%=INKEY(5)
 1250 IF key%=83 phase%=1
 1260 IF key%=69 phase%=3
 1270 UNTIL phase%<>2
 1280 ENDPROC
 1290 :
 1300 DEF PROCplotgalaxies
 1310 PRINT TAB(31,4)a TAB(65,4)d TAB(31,5)b TAB(65,5)e TAB(31,6)c TAB(65,6)f
 1320 PRINT TAB(31,7)o TAB(65,7)r TAB(31,8)p TAB(65,8)s TAB(31,9)q TAB(65,9)t
 1330 PRINT TAB(75,11);T%
 1340 VDU 28,0,31,79,12,12,26
 1350 MOVE 640,0:DRAW 640,640
 1360 MOVE 0,640:DRAW 1279,640
 1370 g=(m*a+n*d)/(m+n)
 1380 h=(m*b+n*e)/(m+n)
 1390 i=(m*c+n*f)/(m+n)
 1400 :
 1410 VDU 26,24,0;0;640;640;
 1420 PROCc((a-g)*4.267+320,(h-b)*4.267+320)
 1430 PROCd((d-g)*4.267+320,(h-e)*4.267+320)
 1440 VDU 26,24,640;0;1279;640;
 1450 PROCc((a-g)*4.267+960,(i-c)*4.267*mz%+320)
 1460 PROCd((d-g)*4.267+960,(i-f)*4.267*mz%+320)
 1470 FOR I%=0 TO stars%-1
 1480 VDU 26,24,0;0;640;640;:PLOT 69,(x(I%)-g)*4.267+320,(h-y(I%))*4.267+320
 1490 VDU 26,24,640;0;1279;640;:PLOT 69,(x(I%)-g)*4.267+960,(i-z(I%))*mz%*4.267+320
 1500 NEXT
 1510 ENDPROC
 1520 :
 1530 DEF PROCc(A%,B%)
 1540 LOCAL I%,a,b
 1550 PLOT 69,A%,B%
 1560 FOR I%=0 TO 3
 1570 a=galc%*COS(PI/8*I%)
 1580 b=galc%*SIN(PI/8*I%)
 1590 PLOT 69,A%+a*2,B%+b*2:PLOT 69,A%-b*2,B%+a*2
 1600 PLOT 69,A%-a*2,B%-b*2:PLOT 69,A%+b*2,B%-a*2
 1610 PLOT 69,A%+a,B%+b:PLOT 69,A%-b,B%+a
 1620 PLOT 69,A%-a,B%-b:PLOT 69,A%+b,B%-a
 1630 NEXT
 1640 ENDPROC
 1650 :
 1660 DEF PROCd(A%,B%)
 1670 LOCAL I%,a
 1680 MOVE A%+galc%,B%
 1690 FOR I%=1 TO 16
 1700 MOVE A%,B%
 1710 a=PI*I%/8
 1720 PLOT 85,A%+galc%*COS(a),B%+galc%*SIN(a)
 1730 NEXT
 1740 ENDPROC
 1750 :
 1760 DEF PROCsdump
 1770 *FX 6,0
 1780 */C.SDUMP
 1790 ENDPROC
 1800 :
 1810 DEF PROCwait
 1820 PRINT TAB(68,10)"..waiting.."
 1830 REPEAT
 1840 key%=INKEY(5)
 1850 UNTIL key%<>-1
 1860 PRINT TAB(68,10)SPC(11)
 1870 ENDPROC
 1880 :
 1890 DEF PROCiterate
 1900 PRINT TAB(11,10)SPC68TAB(11,10)"Iterate until time = ";
 1910 INPUT,stop%
 1920 PROCkey
 1930 ENDPROC
 1940 :
 1950 DEF PROCkey
 1960 COLOUR 129:COLOUR 0
 1970 PRINT TAB(11,10)"(P)ause (I)terations (G)raphics dump (S)et up (E)xit"
 1980 COLOUR 128:COLOUR 1
 1990 ENDPROC
