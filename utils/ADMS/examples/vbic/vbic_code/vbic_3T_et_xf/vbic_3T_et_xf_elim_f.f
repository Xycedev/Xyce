      subroutine vbic_3T_et_xf_f(p,Vrth,Vbei,Vbex,Vbci,Vbep,Vrxf,Vrcx
     +,Vbcx,Vrci,Vrbx,Vrbi,Vre,Vrbp,Vbe,Vbc,Vcei,Vcxf,Ibe,Ibex,Itxf,Itzr
     +,Ibc,Ibep,Ircx,Irci,Irbx,Irbi,Ire,Irbp,Qbe,Qbex,Qbc,Qbcx,Qbep,Qbeo
     +,Qbco,Irth,Ith,Qcth,Ixzf,Ixxf,Qcxf,Flxf,SCALE)
      implicit double precision (a-z)
      double precision p(108)
c
c     Function only code
c
      Tini=2.731500e+02+p(1)
      Tdev=(2.731500e+02+p(1))+Vrth
      Vtv=1.380662e-23*Tdev/1.602189e-19
      rT=Tdev/Tini
      dT=Tdev-Tini
      xvar1=rT**p(91)
      IKFatT=p(54)*xvar1
      xvar1=rT**p(92)
      RCXatT=p(2)*xvar1
      xvar1=rT**p(69)
      RCIatT=p(3)*xvar1
      xvar1=rT**p(93)
      RBXatT=p(7)*xvar1
      xvar1=rT**p(68)
      RBIatT=p(8)*xvar1
      xvar1=rT**p(67)
      REatT=p(9)*xvar1
      xvar1=rT**p(70)
      RSatT=p(10)*xvar1
      xvar1=rT**p(94)
      RBPatT=p(11)*xvar1
      xvar2=rT**p(79)
      xvar3=-p(72)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(13))
      xvar6=xvar1**xvar5
      ISatT=p(12)*xvar6
      xvar2=rT**p(96)
      xvar3=-p(97)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(14))
      xvar6=xvar1**xvar5
      ISRRatT=p(95)*xvar6
      xvar2=rT**p(79)
      xvar3=-p(98)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(45))
      xvar6=xvar1**xvar5
      ISPatT=p(43)*xvar6
      xvar2=rT**p(80)
      xvar3=-p(73)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(34))
      xvar6=xvar1**xvar5
      IBEIatT=p(32)*xvar6
      xvar2=rT**p(81)
      xvar3=-p(76)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(36))
      xvar6=xvar1**xvar5
      IBENatT=p(35)*xvar6
      xvar2=rT**p(80)
      xvar3=-p(74)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(38))
      xvar6=xvar1**xvar5
      IBCIatT=p(37)*xvar6
      xvar2=rT**p(81)
      xvar3=-p(77)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(40))
      xvar6=xvar1**xvar5
      IBCNatT=p(39)*xvar6
      xvar2=rT**p(80)
      xvar3=-p(74)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(38))
      xvar6=xvar1**xvar5
      IBEIPatT=p(46)*xvar6
      xvar2=rT**p(81)
      xvar3=-p(77)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(40))
      xvar6=xvar1**xvar5
      IBENPatT=p(47)*xvar6
      xvar2=rT**p(80)
      xvar3=-p(75)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(49))
      xvar6=xvar1**xvar5
      IBCIPatT=p(48)*xvar6
      xvar2=rT**p(81)
      xvar3=-p(78)*(1.0-rT)/Vtv
      xvar4=exp(xvar3)
      xvar1=(xvar2*xvar4)
      xvar5=(1.0/p(51))
      xvar6=xvar1**xvar5
      IBCNPatT=p(50)*xvar6
      NFatT=p(13)*(1.0+dT*p(82))
      NRatT=p(14)*(1.0+dT*p(82))
      AVC2atT=p(42)*(1.0+dT*p(83))
      VBBEatT=p(99)*(1.0+dT*(p(102)+dT*p(103)))
      NBBEatT=p(100)*(1.0+dT*p(104))
      xvar2=0.5*p(18)*rT/Vtv
      xvar3=exp(xvar2)
      xvar4=-0.5*p(18)*rT/Vtv
      xvar5=exp(xvar4)
      xvar1=xvar3-xvar5
      xvar6=log(xvar1)
      psiio=2.0*(Vtv/rT)*xvar6
      xvar1=log(rT)
      psiin=psiio*rT-3.0*Vtv*xvar1-p(73)*(rT-1.0)
      xvar2=-psiin/Vtv
      xvar3=exp(xvar2)
      xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3))
      xvar4=log(xvar1)
      PEatT=psiin+2.0*Vtv*xvar4
      xvar2=0.5*p(25)*rT/Vtv
      xvar3=exp(xvar2)
      xvar4=-0.5*p(25)*rT/Vtv
      xvar5=exp(xvar4)
      xvar1=xvar3-xvar5
      xvar6=log(xvar1)
      psiio=2.0*(Vtv/rT)*xvar6
      xvar1=log(rT)
      psiin=psiio*rT-3.0*Vtv*xvar1-p(74)*(rT-1.0)
      xvar2=-psiin/Vtv
      xvar3=exp(xvar2)
      xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3))
      xvar4=log(xvar1)
      PCatT=psiin+2.0*Vtv*xvar4
      xvar2=0.5*p(29)*rT/Vtv
      xvar3=exp(xvar2)
      xvar4=-0.5*p(29)*rT/Vtv
      xvar5=exp(xvar4)
      xvar1=xvar3-xvar5
      xvar6=log(xvar1)
      psiio=2.0*(Vtv/rT)*xvar6
      xvar1=log(rT)
      psiin=psiio*rT-3.0*Vtv*xvar1-p(75)*(rT-1.0)
      xvar2=-psiin/Vtv
      xvar3=exp(xvar2)
      xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3))
      xvar4=log(xvar1)
      PSatT=psiin+2.0*Vtv*xvar4
      xvar1=p(18)/PEatT
      xvar2=xvar1**p(19)
      CJEatT=p(17)*xvar2
      xvar1=p(25)/PCatT
      xvar2=xvar1**p(26)
      CJCatT=p(22)*xvar2
      xvar1=p(25)/PCatT
      xvar2=xvar1**p(26)
      CJEPatT=p(24)*xvar2
      xvar1=p(29)/PSatT
      xvar2=xvar1**p(30)
      CJCPatT=p(28)*xvar2
      xvar1=rT**p(79)
      xvar2=-p(72)*(1.0-rT)/Vtv
      xvar3=exp(xvar2)
      GAMMatT=p(5)*xvar1*xvar3
      xvar1=rT**p(71)
      VOatT=p(4)*xvar1
      xvar1=-VBBEatT/(NBBEatT*Vtv)
      EBBEatT=exp(xvar1)
      if(p(52).GT.0.0)then
         IVEF=1.0/p(52)
      else
         IVEF=0.0
      endif
      if(p(53).GT.0.0)then
         IVER=1.0/p(53)
      else
         IVER=0.0
      endif
      if(p(54).GT.0.0)then
         IIKF=1.0/IKFatT
      else
         IIKF=0.0
      endif
      if(p(55).GT.0.0)then
         IIKR=1.0/p(55)
      else
         IIKR=0.0
      endif
      if(p(56).GT.0.0)then
         IIKP=1.0/p(56)
      else
         IIKP=0.0
      endif
      if(p(4).GT.0.0)then
         IVO=1.0/VOatT
      else
         IVO=0.0
      endif
      if(p(6).GT.0.0)then
         IHRCF=1.0/p(6)
      else
         IHRCF=0.0
      endif
      if(p(60).GT.0.0)then
         IVTF=1.0/p(60)
      else
         IVTF=0.0
      endif
      if(p(61).GT.0.0)then
         IITF=1.0/p(61)
      else
         IITF=0.0
      endif
      if(p(61).GT.0.0)then
         slTF=0.0
      else
         slTF=1.0
      endif
      if(p(63).GT.0.0)then
         LEP=p(63)/3.0
      else
         LEP=0.0
      endif
      if(p(63).GT.0.0)then
         CEP=p(63)
      else
         CEP=0.0
      endif
      dv0=-PEatT*p(15)
      if(p(20).LE.0.0)then
         dvh=Vbei+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-p(15))
            xvar2=(-1.0-p(19))
            pwq=xvar1**xvar2
            qlo=PEatT*(1.0-pwq*(1.0-p(15))*(1.0-p(15)))/(1.0-p(19))
            qhi=dvh*(1.0-p(15)+0.5*p(19)*dvh/PEatT)*pwq
         else
            xvar1=(1.0-Vbei/PEatT)
            xvar2=(1.0-p(19))
            xvar3=xvar1**xvar2
            qlo=PEatT*(1.0-xvar3)/(1.0-p(19))
            qhi=0.0
         endif
         qdbe=qlo+qhi
      else
         mv0=sqrt(dv0*dv0+4.0*p(20)*p(20))
         vl0=-0.5*(dv0+mv0)
         xvar1=(1.0-vl0/PEatT)
         xvar2=(1.0-p(19))
         xvar3=xvar1**xvar2
         q0=-PEatT*xvar3/(1.0-p(19))
         dv=Vbei+dv0
         mv=sqrt(dv*dv+4.0*p(20)*p(20))
         vl=0.5*(dv-mv)-dv0
         xvar1=(1.0-vl/PEatT)
         xvar2=(1.0-p(19))
         xvar3=xvar1**xvar2
         qlo=-PEatT*xvar3/(1.0-p(19))
         xvar1=(1.0-p(15))
         xvar2=(-p(19))
         xvar3=xvar1**xvar2
         qdbe=qlo+xvar3*(Vbei-vl+vl0)-q0
      endif
      dv0=-PEatT*p(15)
      if(p(20).LE.0.0)then
         dvh=Vbex+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-p(15))
            xvar2=(-1.0-p(19))
            pwq=xvar1**xvar2
            qlo=PEatT*(1.0-pwq*(1.0-p(15))*(1.0-p(15)))/(1.0-p(19))
            qhi=dvh*(1.0-p(15)+0.5*p(19)*dvh/PEatT)*pwq
         else
            xvar1=(1.0-Vbex/PEatT)
            xvar2=(1.0-p(19))
            xvar3=xvar1**xvar2
            qlo=PEatT*(1.0-xvar3)/(1.0-p(19))
            qhi=0.0
         endif
         qdbex=qlo+qhi
      else
         mv0=sqrt(dv0*dv0+4.0*p(20)*p(20))
         vl0=-0.5*(dv0+mv0)
         xvar1=(1.0-vl0/PEatT)
         xvar2=(1.0-p(19))
         xvar3=xvar1**xvar2
         q0=-PEatT*xvar3/(1.0-p(19))
         dv=Vbex+dv0
         mv=sqrt(dv*dv+4.0*p(20)*p(20))
         vl=0.5*(dv-mv)-dv0
         xvar1=(1.0-vl/PEatT)
         xvar2=(1.0-p(19))
         xvar3=xvar1**xvar2
         qlo=-PEatT*xvar3/(1.0-p(19))
         xvar1=(1.0-p(15))
         xvar2=(-p(19))
         xvar3=xvar1**xvar2
         qdbex=qlo+xvar3*(Vbex-vl+vl0)-q0
      endif
      dv0=-PCatT*p(15)
      if(p(27).LE.0.0)then
         dvh=Vbci+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-p(15))
            xvar2=(-1.0-p(26))
            pwq=xvar1**xvar2
            qlo=PCatT*(1.0-pwq*(1.0-p(15))*(1.0-p(15)))/(1.0-p(26))
            qhi=dvh*(1.0-p(15)+0.5*p(26)*dvh/PCatT)*pwq
         else
            if((p(86).GT.0.0).AND.(Vbci.LT.-p(86)))then
               xvar1=(1.0+p(86)/PCatT)
               xvar2=(1.0-p(26))
               xvar3=xvar1**xvar2
               qlo=PCatT*(1.0-xvar3*(1.0-((1.0-p(26))*(Vbci+p(86)))/(PCa
     +tT+p(86))))/(1.0-p(26))
            else
               xvar1=(1.0-Vbci/PCatT)
               xvar2=(1.0-p(26))
               xvar3=xvar1**xvar2
               qlo=PCatT*(1.0-xvar3)/(1.0-p(26))
            endif
            qhi=0.0
         endif
         qdbc=qlo+qhi
      else
         if((p(86).GT.0.0).AND.(p(87).GT.0.0))then
            vn0=(p(86)+dv0)/(p(86)-dv0)
            vnl0=2.0*vn0/(sqrt((vn0-1.0)*(vn0-1.0)+4.0*p(27)*p(27))+sqrt
     +((vn0+1.0)*(vn0+1.0)+4.0*p(87)*p(87)))
            vl0=0.5*(vnl0*(p(86)-dv0)-p(86)-dv0)
            xvar1=(1.0-vl0/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            qlo0=PCatT*(1.0-xvar3)/(1.0-p(26))
            vn=(2.0*Vbci+p(86)+dv0)/(p(86)-dv0)
            vnl=2.0*vn/(sqrt((vn-1.0)*(vn-1.0)+4.0*p(27)*p(27))+sqrt((vn
     ++1.0)*(vn+1.0)+4.0*p(87)*p(87)))
            vl=0.5*(vnl*(p(86)-dv0)-p(86)-dv0)
            xvar1=(1.0-vl/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            qlo=PCatT*(1.0-xvar3)/(1.0-p(26))
            sel=0.5*(vnl+1.0)
            xvar1=(1.0+p(86)/PCatT)
            xvar2=(-p(26))
            crt=xvar1**xvar2
            xvar1=(1.0+dv0/PCatT)
            xvar2=(-p(26))
            cmx=xvar1**xvar2
            cl=(1.0-sel)*crt+sel*cmx
            ql=(Vbci-vl+vl0)*cl
            qdbc=ql+qlo-qlo0
         else
            mv0=sqrt(dv0*dv0+4.0*p(27)*p(27))
            vl0=-0.5*(dv0+mv0)
            xvar1=(1.0-vl0/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            q0=-PCatT*xvar3/(1.0-p(26))
            dv=Vbci+dv0
            mv=sqrt(dv*dv+4.0*p(27)*p(27))
            vl=0.5*(dv-mv)-dv0
            xvar1=(1.0-vl/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            qlo=-PCatT*xvar3/(1.0-p(26))
            xvar1=(1.0-p(15))
            xvar2=(-p(26))
            xvar3=xvar1**xvar2
            qdbc=qlo+xvar3*(Vbci-vl+vl0)-q0
         endif
      endif
      dv0=-PCatT*p(15)
      if(p(27).LE.0.0)then
         dvh=Vbep+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-p(15))
            xvar2=(-1.0-p(26))
            pwq=xvar1**xvar2
            qlo=PCatT*(1.0-pwq*(1.0-p(15))*(1.0-p(15)))/(1.0-p(26))
            qhi=dvh*(1.0-p(15)+0.5*p(26)*dvh/PCatT)*pwq
         else
            if((p(86).GT.0.0).AND.(Vbep.LT.-p(86)))then
               xvar1=(1.0+p(86)/PCatT)
               xvar2=(1.0-p(26))
               xvar3=xvar1**xvar2
               qlo=PCatT*(1.0-xvar3*(1.0-((1.0-p(26))*(Vbep+p(86)))/(PCa
     +tT+p(86))))/(1.0-p(26))
            else
               xvar1=(1.0-Vbep/PCatT)
               xvar2=(1.0-p(26))
               xvar3=xvar1**xvar2
               qlo=PCatT*(1.0-xvar3)/(1.0-p(26))
            endif
            qhi=0.0
         endif
         qdbep=qlo+qhi
      else
         if((p(86).GT.0.0).AND.(p(87).GT.0.0))then
            vn0=(p(86)+dv0)/(p(86)-dv0)
            vnl0=2.0*vn0/(sqrt((vn0-1.0)*(vn0-1.0)+4.0*p(27)*p(27))+sqrt
     +((vn0+1.0)*(vn0+1.0)+4.0*p(87)*p(87)))
            vl0=0.5*(vnl0*(p(86)-dv0)-p(86)-dv0)
            xvar1=(1.0-vl0/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            qlo0=PCatT*(1.0-xvar3)/(1.0-p(26))
            vn=(2.0*Vbep+p(86)+dv0)/(p(86)-dv0)
            vnl=2.0*vn/(sqrt((vn-1.0)*(vn-1.0)+4.0*p(27)*p(27))+sqrt((vn
     ++1.0)*(vn+1.0)+4.0*p(87)*p(87)))
            vl=0.5*(vnl*(p(86)-dv0)-p(86)-dv0)
            xvar1=(1.0-vl/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            qlo=PCatT*(1.0-xvar3)/(1.0-p(26))
            sel=0.5*(vnl+1.0)
            xvar1=(1.0+p(86)/PCatT)
            xvar2=(-p(26))
            crt=xvar1**xvar2
            xvar1=(1.0+dv0/PCatT)
            xvar2=(-p(26))
            cmx=xvar1**xvar2
            cl=(1.0-sel)*crt+sel*cmx
            ql=(Vbep-vl+vl0)*cl
            qdbep=ql+qlo-qlo0
         else
            mv0=sqrt(dv0*dv0+4.0*p(27)*p(27))
            vl0=-0.5*(dv0+mv0)
            xvar1=(1.0-vl0/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            q0=-PCatT*xvar3/(1.0-p(26))
            dv=Vbep+dv0
            mv=sqrt(dv*dv+4.0*p(27)*p(27))
            vl=0.5*(dv-mv)-dv0
            xvar1=(1.0-vl/PCatT)
            xvar2=(1.0-p(26))
            xvar3=xvar1**xvar2
            qlo=-PCatT*xvar3/(1.0-p(26))
            xvar1=(1.0-p(15))
            xvar2=(-p(26))
            xvar3=xvar1**xvar2
            qdbep=qlo+xvar3*(Vbep-vl+vl0)-q0
         endif
      endif
      argi=Vbei/(NFatT*Vtv)
      if(argi.GT.100.0)then
         expi=exp(100.0)*(1.0+(argi-100.0))
      elseif(argi.LT.-100.0)then
         expi=exp(-100.0)*(1.0+(argi+100.0))
      else
         expi=exp(argi)
      endif
      Ifi=ISatT*(expi-1.0)
      argi=Vbci/(NRatT*Vtv)
      if(argi.GT.100.0)then
         expi=exp(100.0)*(1.0+(argi-100.0))
      elseif(argi.LT.-100.0)then
         expi=exp(-100.0)*(1.0+(argi+100.0))
      else
         expi=exp(argi)
      endif
      Iri=ISatT*ISRRatT*(expi-1.0)
      q1z=1.0+qdbe*IVER+qdbc*IVEF
      q1=0.5*(sqrt((q1z-1.0e-4)*(q1z-1.0e-4)+1.0e-8)+q1z-1.0e-4)+1.0e-4
      q2=Ifi*IIKF+Iri*IIKR
      if(p(89).LT.0.5)then
         xvar2=1.0/p(90)
         xvar3=q1**xvar2
         xvar1=(xvar3+4.0*q2)
         xvar4=xvar1**p(90)
         qb=0.5*(q1+xvar4)
      else
         xvar1=(1.0+4.0*q2)
         xvar2=xvar1**p(90)
         qb=0.5*q1*(1.0+xvar2)
      endif
      Itzr=Iri/qb
      Itzf=Ifi/qb
      Ixzf=-Itzf
      Itxf=Vrxf
      Ixxf=Vrxf
      if(p(43).GT.0.0)then
         argi=Vbep/(p(45)*Vtv)
         if(argi.GT.100.0)then
            expi=exp(100.0)*(1.0+(argi-100.0))
         elseif(argi.LT.-100.0)then
            expi=exp(-100.0)*(1.0+(argi+100.0))
         else
            expi=exp(argi)
         endif
         argx=Vbci/(p(45)*Vtv)
         if(argx.GT.100.0)then
            expx=exp(100.0)*(1.0+(argx-100.0))
         elseif(argx.LT.-100.0)then
            expx=exp(-100.0)*(1.0+(argx+100.0))
         else
            expx=exp(argx)
         endif
         Ifp=ISPatT*(p(44)*expi+(1.0-p(44))*expx-1.0)
         q2p=Ifp*IIKP
         qbp=0.5*(1.0+sqrt(1.0+4.0*q2p))
      else
         Ifp=0.0
         qbp=1.0
      endif
      if(p(33).EQ.1.0)then
         argi=Vbei/(p(34)*Vtv)
         if(argi.GT.100.0)then
            expi=exp(100.0)*(1.0+(argi-100.0))
         elseif(argi.LT.-100.0)then
            expi=exp(-100.0)*(1.0+(argi+100.0))
         else
            expi=exp(argi)
         endif
         argn=Vbei/(p(36)*Vtv)
         if(argn.GT.100.0)then
            expn=exp(100.0)*(1.0+(argn-100.0))
         elseif(argn.LT.-100.0)then
            expn=exp(-100.0)*(1.0+(argn+100.0))
         else
            expn=exp(argn)
         endif
         if(p(99).GT.0.0)then
            argx=(-VBBEatT-Vbei)/(NBBEatT*Vtv)
            if(argx.GT.100.0)then
               expx=exp(100.0)*(1.0+(argx-100.0))
            elseif(argx.LT.-100.0)then
               expx=exp(-100.0)*(1.0+(argx+100.0))
            else
               expx=exp(argx)
            endif
            Ibe=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p(101)*(expx-EBBEa
     +tT)
         else
            Ibe=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)
         endif
         Ibex=0.0
      elseif(p(33).EQ.0.0)then
         Ibe=0.0
         argi=Vbex/(p(34)*Vtv)
         if(argi.GT.100.0)then
            expi=exp(100.0)*(1.0+(argi-100.0))
         elseif(argi.LT.-100.0)then
            expi=exp(-100.0)*(1.0+(argi+100.0))
         else
            expi=exp(argi)
         endif
         argn=Vbex/(p(36)*Vtv)
         if(argn.GT.100.0)then
            expn=exp(100.0)*(1.0+(argn-100.0))
         elseif(argn.LT.-100.0)then
            expn=exp(-100.0)*(1.0+(argn+100.0))
         else
            expn=exp(argn)
         endif
         if(p(99).GT.0.0)then
            argx=(-VBBEatT-Vbex)/(NBBEatT*Vtv)
            if(argx.GT.100.0)then
               expx=exp(100.0)*(1.0+(argx-100.0))
            elseif(argx.LT.-100.0)then
               expx=exp(-100.0)*(1.0+(argx+100.0))
            else
               expx=exp(argx)
            endif
            Ibex=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p(101)*(expx-EBBE
     +atT)
         else
            Ibex=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)
         endif
      else
         argi=Vbei/(p(34)*Vtv)
         if(argi.GT.100.0)then
            expi=exp(100.0)*(1.0+(argi-100.0))
         elseif(argi.LT.-100.0)then
            expi=exp(-100.0)*(1.0+(argi+100.0))
         else
            expi=exp(argi)
         endif
         argn=Vbei/(p(36)*Vtv)
         if(argn.GT.100.0)then
            expn=exp(100.0)*(1.0+(argn-100.0))
         elseif(argn.LT.-100.0)then
            expn=exp(-100.0)*(1.0+(argn+100.0))
         else
            expn=exp(argn)
         endif
         if(p(99).GT.0.0)then
            argx=(-VBBEatT-Vbei)/(NBBEatT*Vtv)
            if(argx.GT.100.0)then
               expx=exp(100.0)*(1.0+(argx-100.0))
            elseif(argx.LT.-100.0)then
               expx=exp(-100.0)*(1.0+(argx+100.0))
            else
               expx=exp(argx)
            endif
            Ibe=p(33)*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p(101)*(exp
     +x-EBBEatT))
         else
            Ibe=p(33)*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0))
         endif
         argi=Vbex/(p(34)*Vtv)
         if(argi.GT.100.0)then
            expi=exp(100.0)*(1.0+(argi-100.0))
         elseif(argi.LT.-100.0)then
            expi=exp(-100.0)*(1.0+(argi+100.0))
         else
            expi=exp(argi)
         endif
         argn=Vbex/(p(36)*Vtv)
         if(argn.GT.100.0)then
            expn=exp(100.0)*(1.0+(argn-100.0))
         elseif(argn.LT.-100.0)then
            expn=exp(-100.0)*(1.0+(argn+100.0))
         else
            expn=exp(argn)
         endif
         if(p(99).GT.0.0)then
            argx=(-VBBEatT-Vbex)/(NBBEatT*Vtv)
            if(argx.GT.100.0)then
               expx=exp(100.0)*(1.0+(argx-100.0))
            elseif(argx.LT.-100.0)then
               expx=exp(-100.0)*(1.0+(argx+100.0))
            else
               expx=exp(argx)
            endif
            Ibex=(1.0-p(33))*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p(10
     +1)*(expx-EBBEatT))
         else
            Ibex=(1.0-p(33))*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0))
         endif
      endif
      argi=Vbci/(p(38)*Vtv)
      if(argi.GT.100.0)then
         expi=exp(100.0)*(1.0+(argi-100.0))
      elseif(argi.LT.-100.0)then
         expi=exp(-100.0)*(1.0+(argi+100.0))
      else
         expi=exp(argi)
      endif
      argn=Vbci/(p(40)*Vtv)
      if(argn.GT.100.0)then
         expn=exp(100.0)*(1.0+(argn-100.0))
      elseif(argn.LT.-100.0)then
         expn=exp(-100.0)*(1.0+(argn+100.0))
      else
         expn=exp(argn)
      endif
      Ibcj=IBCIatT*(expi-1.0)+IBCNatT*(expn-1.0)
      if((p(46).GT.0.0).OR.(p(47).GT.0.0))then
         argi=Vbep/(p(38)*Vtv)
         if(argi.GT.100.0)then
            expi=exp(100.0)*(1.0+(argi-100.0))
         elseif(argi.LT.-100.0)then
            expi=exp(-100.0)*(1.0+(argi+100.0))
         else
            expi=exp(argi)
         endif
         argn=Vbep/(p(40)*Vtv)
         if(argn.GT.100.0)then
            expn=exp(100.0)*(1.0+(argn-100.0))
         elseif(argn.LT.-100.0)then
            expn=exp(-100.0)*(1.0+(argn+100.0))
         else
            expn=exp(argn)
         endif
         Ibep=IBEIPatT*(expi-1.0)+IBENPatT*(expn-1.0)
      else
         Ibep=0.0
      endif
      if(p(41).GT.0.0)then
         vl=0.5*(sqrt((PCatT-Vbci)*(PCatT-Vbci)+0.01)+(PCatT-Vbci))
         xvar2=(p(26)-1.0)
         xvar3=vl**xvar2
         xvar1=-AVC2atT*xvar3
         if(xvar1.GT.100.0)then
            xvar4=exp(100.0)*(1.0+(xvar1-100.0))
         elseif(xvar1.LT.-100.0)then
            xvar4=exp(-100.0)*(1.0+(xvar1+100.0))
         else
            xvar4=exp(xvar1)
         endif
         avalf=p(41)*vl*xvar4
         Igc=(Itxf-Itzr-Ibcj)*avalf
      else
         Igc=0.0
      endif
      Ibc=Ibcj-Igc
      if(p(2).GT.0.0)then
         Ircx=Vrcx/RCXatT
      else
         Ircx=0.0
      endif
      argi=Vbci/Vtv
      if(argi.GT.100.0)then
         expi=exp(100.0)*(1.0+(argi-100.0))
      elseif(argi.LT.-100.0)then
         expi=exp(-100.0)*(1.0+(argi+100.0))
      else
         expi=exp(argi)
      endif
      argx=Vbcx/Vtv
      if(argx.GT.100.0)then
         expx=exp(100.0)*(1.0+(argx-100.0))
      elseif(argx.LT.-100.0)then
         expx=exp(-100.0)*(1.0+(argx+100.0))
      else
         expx=exp(argx)
      endif
      Kbci=sqrt(1.0+GAMMatT*expi)
      Kbcx=sqrt(1.0+GAMMatT*expx)
      if(p(3).GT.0.0)then
         rKp1=(Kbci+1.0)/(Kbcx+1.0)
         xvar1=log(rKp1)
         Iohm=(Vrci+Vtv*(Kbci-Kbcx-xvar1))/RCIatT
         derf=IVO*RCIatT*Iohm/(1.0+0.5*IVO*IHRCF*sqrt(Vrci*Vrci+0.01))
         Irci=Iohm/sqrt(1.0+derf*derf)
      else
         Irci=0.0
      endif
      if(p(7).GT.0.0)then
         Irbx=Vrbx/RBXatT
      else
         Irbx=0.0
      endif
      if(p(8).GT.0.0)then
         Irbi=Vrbi*qb/RBIatT
      else
         Irbi=0.0
      endif
      if(p(9).GT.0.0)then
         Ire=Vre/REatT
      else
         Ire=0.0
      endif
      if(p(11).GT.0.0)then
         Irbp=Vrbp*qbp/RBPatT
      else
         Irbp=0.0
      endif
      if(Ifi.GT.0.0)then
         sgIf=1.0
      else
         sgIf=0.0
      endif
      rIf=Ifi*sgIf*IITF
      mIf=rIf/(rIf+1.0)
      xvar1=Vbci*IVTF/1.44
      if(xvar1.GT.100.0)then
         xvar2=exp(100.0)*(1.0+(xvar1-100.0))
      elseif(xvar1.LT.-100.0)then
         xvar2=exp(-100.0)*(1.0+(xvar1+100.0))
      else
         xvar2=exp(xvar1)
      endif
      tff=p(57)*(1.0+p(58)*q1)*(1.0+p(59)*xvar2*(slTF+mIf*mIf)*sgIf)
      Qbe=CJEatT*qdbe*p(33)+tff*Ifi/qb
      Qbex=CJEatT*qdbex*(1.0-p(33))
      Qbc=CJCatT*qdbc+p(62)*Iri+p(23)*Kbci
      Qbcx=p(23)*Kbcx
      Qbep=CJEPatT*qdbep+p(62)*Ifp
      Qbeo=Vbe*p(16)
      Qbco=Vbc*p(21)
      Ith=-(Ibe*Vbei+Ibc*Vbci+(Itxf-Itzr)*Vcei+Ibex*Vbex+Ibep*Vbep+Ircx*
     +Vrcx+Irci*Vrci+Irbx*Vrbx+Irbi*Vrbi+Ire*Vre+Irbp*Vrbp)
      if(p(84).GT.0.0)then
         Irth=Vrth/p(84)
      else
         Irth=0.0
      endif
      Qcth=Vrth*p(85)
      Flxf=LEP*Ixxf
      Qcxf=CEP*Vcxf
c
c     Scale outputs
c
      if(SCALE.NE.1.0)then
         Ibe=SCALE*Ibe
         Ibex=SCALE*Ibex
         Itxf=SCALE*Itxf
         Itzr=SCALE*Itzr
         Ibc=SCALE*Ibc
         Ibep=SCALE*Ibep
         Ircx=SCALE*Ircx
         Irci=SCALE*Irci
         Irbx=SCALE*Irbx
         Irbi=SCALE*Irbi
         Ire=SCALE*Ire
         Irbp=SCALE*Irbp
         Qbe=SCALE*Qbe
         Qbex=SCALE*Qbex
         Qbc=SCALE*Qbc
         Qbcx=SCALE*Qbcx
         Qbep=SCALE*Qbep
         Qbeo=SCALE*Qbeo
         Qbco=SCALE*Qbco
         Irth=SCALE*Irth
         Ith=SCALE*Ith
         Qcth=SCALE*Qcth
         Ixzf=SCALE*Ixzf
         Ixxf=SCALE*Ixxf
         Qcxf=SCALE*Qcxf
         Flxf=SCALE*Flxf
      endif
      return
      end
