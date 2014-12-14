      subroutine vbic_4T_it_cf_t(p,pnom,TAMB)
      implicit double precision (a-z)
      double precision p(108),pnom(108)
c
c     Direct copy p<-pnom for temperature independent parameters
c
      p(6)=pnom(6)
      p(15)=pnom(15)
      p(16)=pnom(16)
      p(19)=pnom(19)
      p(20)=pnom(20)
      p(21)=pnom(21)
      p(23)=pnom(23)
      p(26)=pnom(26)
      p(27)=pnom(27)
      p(30)=pnom(30)
      p(31)=pnom(31)
      p(33)=pnom(33)
      p(34)=pnom(34)
      p(36)=pnom(36)
      p(38)=pnom(38)
      p(40)=pnom(40)
      p(41)=pnom(41)
      p(44)=pnom(44)
      p(45)=pnom(45)
      p(49)=pnom(49)
      p(51)=pnom(51)
      p(52)=pnom(52)
      p(53)=pnom(53)
      p(55)=pnom(55)
      p(56)=pnom(56)
      p(57)=pnom(57)
      p(58)=pnom(58)
      p(59)=pnom(59)
      p(60)=pnom(60)
      p(61)=pnom(61)
      p(62)=pnom(62)
      p(63)=pnom(63)
      p(64)=pnom(64)
      p(65)=pnom(65)
      p(66)=pnom(66)
      p(67)=pnom(67)
      p(68)=pnom(68)
      p(69)=pnom(69)
      p(70)=pnom(70)
      p(71)=pnom(71)
      p(72)=pnom(72)
      p(73)=pnom(73)
      p(74)=pnom(74)
      p(75)=pnom(75)
      p(76)=pnom(76)
      p(77)=pnom(77)
      p(78)=pnom(78)
      p(79)=pnom(79)
      p(80)=pnom(80)
      p(81)=pnom(81)
      p(82)=pnom(82)
      p(83)=pnom(83)
      p(84)=pnom(84)
      p(85)=pnom(85)
      p(86)=pnom(86)
      p(87)=pnom(87)
      p(88)=pnom(88)
      p(89)=pnom(89)
      p(90)=pnom(90)
      p(91)=pnom(91)
      p(92)=pnom(92)
      p(93)=pnom(93)
      p(94)=pnom(94)
      p(96)=pnom(96)
      p(97)=pnom(97)
      p(98)=pnom(98)
      p(101)=pnom(101)
      p(102)=pnom(102)
      p(103)=pnom(103)
      p(104)=pnom(104)
      p(106)=pnom(106)
      p(107)=pnom(107)
      p(108)=pnom(108)
c
c     Temperature mappings for model parameters
c
      Tini=2.731500e+02+pnom(1)
      Tdev=(2.731500e+02+TAMB)+pnom(106)
      Vtv=1.380662e-23*Tdev/1.602189e-19
      rT=Tdev/Tini
      dT=Tdev-Tini
      xvar1=rT**pnom(91)
      p(54)=pnom(54)*xvar1
      xvar1=rT**pnom(92)
      p(2)=pnom(2)*xvar1
      xvar1=rT**pnom(69)
      p(3)=pnom(3)*xvar1
      xvar1=rT**pnom(93)
      p(7)=pnom(7)*xvar1
      xvar1=rT**pnom(68)
      p(8)=pnom(8)*xvar1
      xvar1=rT**pnom(67)
      p(9)=pnom(9)*xvar1
      xvar1=rT**pnom(70)
      p(10)=pnom(10)*xvar1
      xvar1=rT**pnom(94)
      p(11)=pnom(11)*xvar1
      xvar2=rT**pnom(79)
      xvar3=-pnom(72)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(13))
      xvar6=xvar1**xvar5
      p(12)=pnom(12)*xvar6
      xvar2=rT**pnom(96)
      xvar3=-pnom(97)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(14))
      xvar6=xvar1**xvar5
      p(95)=pnom(95)*xvar6
      xvar2=rT**pnom(79)
      xvar3=-pnom(98)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(45))
      xvar6=xvar1**xvar5
      p(43)=pnom(43)*xvar6
      xvar2=rT**pnom(80)
      xvar3=-pnom(73)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(34))
      xvar6=xvar1**xvar5
      p(32)=pnom(32)*xvar6
      xvar2=rT**pnom(81)
      xvar3=-pnom(76)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(36))
      xvar6=xvar1**xvar5
      p(35)=pnom(35)*xvar6
      xvar2=rT**pnom(80)
      xvar3=-pnom(74)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(38))
      xvar6=xvar1**xvar5
      p(37)=pnom(37)*xvar6
      xvar2=rT**pnom(81)
      xvar3=-pnom(77)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(40))
      xvar6=xvar1**xvar5
      p(39)=pnom(39)*xvar6
      xvar2=rT**pnom(80)
      xvar3=-pnom(74)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(38))
      xvar6=xvar1**xvar5
      p(46)=pnom(46)*xvar6
      xvar2=rT**pnom(81)
      xvar3=-pnom(77)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(40))
      xvar6=xvar1**xvar5
      p(47)=pnom(47)*xvar6
      xvar2=rT**pnom(80)
      xvar3=-pnom(75)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(49))
      xvar6=xvar1**xvar5
      p(48)=pnom(48)*xvar6
      xvar2=rT**pnom(81)
      xvar3=-pnom(78)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/pnom(51))
      xvar6=xvar1**xvar5
      p(50)=pnom(50)*xvar6
      p(13)=pnom(13)*(1.0+dT*pnom(82))
      p(14)=pnom(14)*(1.0+dT*pnom(82))
      p(42)=pnom(42)*(1.0+dT*pnom(83))
      p(99)=pnom(99)*(1.0+dT*(pnom(102)+dT*pnom(103)))
      p(100)=pnom(100)*(1.0+dT*pnom(104))
      xvar2=0.5*pnom(18)*rT/Vtv
      xvar3=exp(xvar2)
      xvar4=-0.5*pnom(18)*rT/Vtv
      xvar5=exp(xvar4)
      xvar1=xvar3-xvar5
      xvar6=log(xvar1)
      psiio=2.0*(Vtv/rT)*xvar6
      xvar1=log(rT)
      psiin=psiio*rT-3.0*Vtv*xvar1-pnom(73)*(rT-1.0)
      xvar2=-psiin/Vtv
      xvar3=exp(xvar2)
      xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3))
      xvar4=log(xvar1)
      p(18)=psiin+2.0*Vtv*xvar4
      xvar2=0.5*pnom(25)*rT/Vtv
      xvar3=exp(xvar2)
      xvar4=-0.5*pnom(25)*rT/Vtv
      xvar5=exp(xvar4)
      xvar1=xvar3-xvar5
      xvar6=log(xvar1)
      psiio=2.0*(Vtv/rT)*xvar6
      xvar1=log(rT)
      psiin=psiio*rT-3.0*Vtv*xvar1-pnom(74)*(rT-1.0)
      xvar2=-psiin/Vtv
      xvar3=exp(xvar2)
      xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3))
      xvar4=log(xvar1)
      p(25)=psiin+2.0*Vtv*xvar4
      xvar2=0.5*pnom(29)*rT/Vtv
      xvar3=exp(xvar2)
      xvar4=-0.5*pnom(29)*rT/Vtv
      xvar5=exp(xvar4)
      xvar1=xvar3-xvar5
      xvar6=log(xvar1)
      psiio=2.0*(Vtv/rT)*xvar6
      xvar1=log(rT)
      psiin=psiio*rT-3.0*Vtv*xvar1-pnom(75)*(rT-1.0)
      xvar2=-psiin/Vtv
      xvar3=exp(xvar2)
      xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3))
      xvar4=log(xvar1)
      p(29)=psiin+2.0*Vtv*xvar4
      xvar1=pnom(18)/p(18)
      xvar2=xvar1**pnom(19)
      p(17)=pnom(17)*xvar2
      xvar1=pnom(25)/p(25)
      xvar2=xvar1**pnom(26)
      p(22)=pnom(22)*xvar2
      xvar1=pnom(25)/p(25)
      xvar2=xvar1**pnom(26)
      p(24)=pnom(24)*xvar2
      xvar1=pnom(29)/p(29)
      xvar2=xvar1**pnom(30)
      p(28)=pnom(28)*xvar2
      xvar1=rT**pnom(79)
      xvar2=-pnom(72)*(1.0-rT)/Vtv
      xvar3=exp(xvar2)
      p(5)=pnom(5)*xvar1*xvar3
      xvar1=rT**pnom(71)
      p(4)=pnom(4)*xvar1
      xvar1=-p(99)/(p(100)*Vtv)
      p(105)=exp(xvar1)
      p(1)=TAMB+p(106)
      return
      end
