      subroutine vbic_4T_it_cf_f(p,Vbei,Vbex,Vbci,Vbep,Vbcp,Vrcx,Vbcx
     +,Vrci,Vrbx,Vrbi,Vre,Vrbp,Vrs,Vbe,Vbc,Ibe,Ibex,Itzf,Itzr,Ibc,Ibep
     +,Ircx,Irci,Irbx,Irbi,Ire,Irbp,Qbe,Qbex,Qbc,Qbcx,Qbep,Qbeo,Qbco
     +,Ibcp,Iccp,Irs,Qbcp,SCALE)
      implicit double precision (a-z)
      double precision p(108)
c
c     Map parameter array into symbolic names
c
      TNOM=p(1)
      RCX=p(2)
      RCI=p(3)
      VO=p(4)
      GAMM=p(5)
      HRCF=p(6)
      RBX=p(7)
      RBI=p(8)
      RE=p(9)
      RS=p(10)
      RBP=p(11)
      IS=p(12)
      NF=p(13)
      NR=p(14)
      FC=p(15)
      CBEO=p(16)
      CJE=p(17)
      PE=p(18)
      ME=p(19)
      AJE=p(20)
      CBCO=p(21)
      CJC=p(22)
      QCO=p(23)
      CJEP=p(24)
      PC=p(25)
      MC=p(26)
      AJC=p(27)
      CJCP=p(28)
      PS=p(29)
      MS=p(30)
      AJS=p(31)
      IBEI=p(32)
      WBE=p(33)
      NEI=p(34)
      IBEN=p(35)
      NEN=p(36)
      IBCI=p(37)
      NCI=p(38)
      IBCN=p(39)
      NCN=p(40)
      AVC1=p(41)
      AVC2=p(42)
      ISP=p(43)
      WSP=p(44)
      NFP=p(45)
      IBEIP=p(46)
      IBENP=p(47)
      IBCIP=p(48)
      NCIP=p(49)
      IBCNP=p(50)
      NCNP=p(51)
      VEF=p(52)
      VER=p(53)
      IKF=p(54)
      IKR=p(55)
      IKP=p(56)
      TF=p(57)
      QTF=p(58)
      XTF=p(59)
      VTF=p(60)
      ITF=p(61)
      TR=p(62)
      VRT=p(86)
      ART=p(87)
      CCSO=p(88)
      QBM=p(89)
      NKF=p(90)
      ISRR=p(95)
      VBBE=p(99)
      NBBE=p(100)
      IBBE=p(101)
      EBBE=p(105)
c
c     Function only code
c
      Vtv=1.380662e-23*(2.731500e+02+TNOM)/1.602189e-19
      if(VEF.GT.0.0)then
         IVEF=1.0/VEF
      else
         IVEF=0.0
      endif
      if(VER.GT.0.0)then
         IVER=1.0/VER
      else
         IVER=0.0
      endif
      if(IKF.GT.0.0)then
         IIKF=1.0/IKF
      else
         IIKF=0.0
      endif
      if(IKR.GT.0.0)then
         IIKR=1.0/IKR
      else
         IIKR=0.0
      endif
      if(IKP.GT.0.0)then
         IIKP=1.0/IKP
      else
         IIKP=0.0
      endif
      if(VO.GT.0.0)then
         IVO=1.0/VO
      else
         IVO=0.0
      endif
      if(HRCF.GT.0.0)then
         IHRCF=1.0/HRCF
      else
         IHRCF=0.0
      endif
      if(VTF.GT.0.0)then
         IVTF=1.0/VTF
      else
         IVTF=0.0
      endif
      if(ITF.GT.0.0)then
         IITF=1.0/ITF
      else
         IITF=0.0
      endif
      if(ITF.GT.0.0)then
         slTF=0.0
      else
         slTF=1.0
      endif
      dv0=-PE*FC
      if(AJE.LE.0.0)then
         dvh=Vbei+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-FC)
            xvar2=(-1.0-ME)
            pwq=xvar1**xvar2
            qlo=PE*(1.0-pwq*(1.0-FC)*(1.0-FC))/(1.0-ME)
            qhi=dvh*(1.0-FC+0.5*ME*dvh/PE)*pwq
         else
            xvar1=(1.0-Vbei/PE)
            xvar2=(1.0-ME)
            xvar3=xvar1**xvar2
            qlo=PE*(1.0-xvar3)/(1.0-ME)
            qhi=0.0
         endif
         qdbe=qlo+qhi
      else
         mv0=sqrt(dv0*dv0+4.0*AJE*AJE)
         vl0=-0.5*(dv0+mv0)
         xvar1=(1.0-vl0/PE)
         xvar2=(1.0-ME)
         xvar3=xvar1**xvar2
         q0=-PE*xvar3/(1.0-ME)
         dv=Vbei+dv0
         mv=sqrt(dv*dv+4.0*AJE*AJE)
         vl=0.5*(dv-mv)-dv0
         xvar1=(1.0-vl/PE)
         xvar2=(1.0-ME)
         xvar3=xvar1**xvar2
         qlo=-PE*xvar3/(1.0-ME)
         xvar1=(1.0-FC)
         xvar2=(-ME)
         xvar3=xvar1**xvar2
         qdbe=qlo+xvar3*(Vbei-vl+vl0)-q0
      endif
      dv0=-PE*FC
      if(AJE.LE.0.0)then
         dvh=Vbex+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-FC)
            xvar2=(-1.0-ME)
            pwq=xvar1**xvar2
            qlo=PE*(1.0-pwq*(1.0-FC)*(1.0-FC))/(1.0-ME)
            qhi=dvh*(1.0-FC+0.5*ME*dvh/PE)*pwq
         else
            xvar1=(1.0-Vbex/PE)
            xvar2=(1.0-ME)
            xvar3=xvar1**xvar2
            qlo=PE*(1.0-xvar3)/(1.0-ME)
            qhi=0.0
         endif
         qdbex=qlo+qhi
      else
         mv0=sqrt(dv0*dv0+4.0*AJE*AJE)
         vl0=-0.5*(dv0+mv0)
         xvar1=(1.0-vl0/PE)
         xvar2=(1.0-ME)
         xvar3=xvar1**xvar2
         q0=-PE*xvar3/(1.0-ME)
         dv=Vbex+dv0
         mv=sqrt(dv*dv+4.0*AJE*AJE)
         vl=0.5*(dv-mv)-dv0
         xvar1=(1.0-vl/PE)
         xvar2=(1.0-ME)
         xvar3=xvar1**xvar2
         qlo=-PE*xvar3/(1.0-ME)
         xvar1=(1.0-FC)
         xvar2=(-ME)
         xvar3=xvar1**xvar2
         qdbex=qlo+xvar3*(Vbex-vl+vl0)-q0
      endif
      dv0=-PC*FC
      if(AJC.LE.0.0)then
         dvh=Vbci+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-FC)
            xvar2=(-1.0-MC)
            pwq=xvar1**xvar2
            qlo=PC*(1.0-pwq*(1.0-FC)*(1.0-FC))/(1.0-MC)
            qhi=dvh*(1.0-FC+0.5*MC*dvh/PC)*pwq
         else
            if((VRT.GT.0.0).AND.(Vbci.LT.-VRT))then
               xvar1=(1.0+VRT/PC)
               xvar2=(1.0-MC)
               xvar3=xvar1**xvar2
               qlo=PC*(1.0-xvar3*(1.0-((1.0-MC)*(Vbci+VRT))/(PC+VRT)))/(
     +1.0-MC)
            else
               xvar1=(1.0-Vbci/PC)
               xvar2=(1.0-MC)
               xvar3=xvar1**xvar2
               qlo=PC*(1.0-xvar3)/(1.0-MC)
            endif
            qhi=0.0
         endif
         qdbc=qlo+qhi
      else
         if((VRT.GT.0.0).AND.(ART.GT.0.0))then
            vn0=(VRT+dv0)/(VRT-dv0)
            vnl0=2.0*vn0/(sqrt((vn0-1.0)*(vn0-1.0)+4.0*AJC*AJC)+sqrt((vn
     +0+1.0)*(vn0+1.0)+4.0*ART*ART))
            vl0=0.5*(vnl0*(VRT-dv0)-VRT-dv0)
            xvar1=(1.0-vl0/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            qlo0=PC*(1.0-xvar3)/(1.0-MC)
            vn=(2.0*Vbci+VRT+dv0)/(VRT-dv0)
            vnl=2.0*vn/(sqrt((vn-1.0)*(vn-1.0)+4.0*AJC*AJC)+sqrt((vn+1.0
     +)*(vn+1.0)+4.0*ART*ART))
            vl=0.5*(vnl*(VRT-dv0)-VRT-dv0)
            xvar1=(1.0-vl/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            qlo=PC*(1.0-xvar3)/(1.0-MC)
            sel=0.5*(vnl+1.0)
            xvar1=(1.0+VRT/PC)
            xvar2=(-MC)
            crt=xvar1**xvar2
            xvar1=(1.0+dv0/PC)
            xvar2=(-MC)
            cmx=xvar1**xvar2
            cl=(1.0-sel)*crt+sel*cmx
            ql=(Vbci-vl+vl0)*cl
            qdbc=ql+qlo-qlo0
         else
            mv0=sqrt(dv0*dv0+4.0*AJC*AJC)
            vl0=-0.5*(dv0+mv0)
            xvar1=(1.0-vl0/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            q0=-PC*xvar3/(1.0-MC)
            dv=Vbci+dv0
            mv=sqrt(dv*dv+4.0*AJC*AJC)
            vl=0.5*(dv-mv)-dv0
            xvar1=(1.0-vl/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            qlo=-PC*xvar3/(1.0-MC)
            xvar1=(1.0-FC)
            xvar2=(-MC)
            xvar3=xvar1**xvar2
            qdbc=qlo+xvar3*(Vbci-vl+vl0)-q0
         endif
      endif
      dv0=-PC*FC
      if(AJC.LE.0.0)then
         dvh=Vbep+dv0
         if(dvh.GT.0.0)then
            xvar1=(1.0-FC)
            xvar2=(-1.0-MC)
            pwq=xvar1**xvar2
            qlo=PC*(1.0-pwq*(1.0-FC)*(1.0-FC))/(1.0-MC)
            qhi=dvh*(1.0-FC+0.5*MC*dvh/PC)*pwq
         else
            if((VRT.GT.0.0).AND.(Vbep.LT.-VRT))then
               xvar1=(1.0+VRT/PC)
               xvar2=(1.0-MC)
               xvar3=xvar1**xvar2
               qlo=PC*(1.0-xvar3*(1.0-((1.0-MC)*(Vbep+VRT))/(PC+VRT)))/(
     +1.0-MC)
            else
               xvar1=(1.0-Vbep/PC)
               xvar2=(1.0-MC)
               xvar3=xvar1**xvar2
               qlo=PC*(1.0-xvar3)/(1.0-MC)
            endif
            qhi=0.0
         endif
         qdbep=qlo+qhi
      else
         if((VRT.GT.0.0).AND.(ART.GT.0.0))then
            vn0=(VRT+dv0)/(VRT-dv0)
            vnl0=2.0*vn0/(sqrt((vn0-1.0)*(vn0-1.0)+4.0*AJC*AJC)+sqrt((vn
     +0+1.0)*(vn0+1.0)+4.0*ART*ART))
            vl0=0.5*(vnl0*(VRT-dv0)-VRT-dv0)
            xvar1=(1.0-vl0/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            qlo0=PC*(1.0-xvar3)/(1.0-MC)
            vn=(2.0*Vbep+VRT+dv0)/(VRT-dv0)
            vnl=2.0*vn/(sqrt((vn-1.0)*(vn-1.0)+4.0*AJC*AJC)+sqrt((vn+1.0
     +)*(vn+1.0)+4.0*ART*ART))
            vl=0.5*(vnl*(VRT-dv0)-VRT-dv0)
            xvar1=(1.0-vl/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            qlo=PC*(1.0-xvar3)/(1.0-MC)
            sel=0.5*(vnl+1.0)
            xvar1=(1.0+VRT/PC)
            xvar2=(-MC)
            crt=xvar1**xvar2
            xvar1=(1.0+dv0/PC)
            xvar2=(-MC)
            cmx=xvar1**xvar2
            cl=(1.0-sel)*crt+sel*cmx
            ql=(Vbep-vl+vl0)*cl
            qdbep=ql+qlo-qlo0
         else
            mv0=sqrt(dv0*dv0+4.0*AJC*AJC)
            vl0=-0.5*(dv0+mv0)
            xvar1=(1.0-vl0/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            q0=-PC*xvar3/(1.0-MC)
            dv=Vbep+dv0
            mv=sqrt(dv*dv+4.0*AJC*AJC)
            vl=0.5*(dv-mv)-dv0
            xvar1=(1.0-vl/PC)
            xvar2=(1.0-MC)
            xvar3=xvar1**xvar2
            qlo=-PC*xvar3/(1.0-MC)
            xvar1=(1.0-FC)
            xvar2=(-MC)
            xvar3=xvar1**xvar2
            qdbep=qlo+xvar3*(Vbep-vl+vl0)-q0
         endif
      endif
      if(CJCP.GT.0.0)then
         dv0=-PS*FC
         if(AJS.LE.0.0)then
            dvh=Vbcp+dv0
            if(dvh.GT.0.0)then
               xvar1=(1.0-FC)
               xvar2=(-1.0-MS)
               pwq=xvar1**xvar2
               qlo=PS*(1.0-pwq*(1.0-FC)*(1.0-FC))/(1.0-MS)
               qhi=dvh*(1.0-FC+0.5*MS*dvh/PS)*pwq
            else
               xvar1=(1.0-Vbcp/PS)
               xvar2=(1.0-MS)
               xvar3=xvar1**xvar2
               qlo=PS*(1.0-xvar3)/(1.0-MS)
               qhi=0.0
            endif
            qdbcp=qlo+qhi
         else
            mv0=sqrt(dv0*dv0+4.0*AJS*AJS)
            vl0=-0.5*(dv0+mv0)
            xvar1=(1.0-vl0/PS)
            xvar2=(1.0-MS)
            xvar3=xvar1**xvar2
            q0=-PS*xvar3/(1.0-MS)
            dv=Vbcp+dv0
            mv=sqrt(dv*dv+4.0*AJS*AJS)
            vl=0.5*(dv-mv)-dv0
            xvar1=(1.0-vl/PS)
            xvar2=(1.0-MS)
            xvar3=xvar1**xvar2
            qlo=-PS*xvar3/(1.0-MS)
            xvar1=(1.0-FC)
            xvar2=(-MS)
            xvar3=xvar1**xvar2
            qdbcp=qlo+xvar3*(Vbcp-vl+vl0)-q0
         endif
      else
         qdbcp=0.0
      endif
      argi=Vbei/(NF*Vtv)
      expi=exp(argi)
      Ifi=IS*(expi-1.0)
      argi=Vbci/(NR*Vtv)
      expi=exp(argi)
      Iri=IS*ISRR*(expi-1.0)
      q1z=1.0+qdbe*IVER+qdbc*IVEF
      q1=0.5*(sqrt((q1z-1.0e-4)*(q1z-1.0e-4)+1.0e-8)+q1z-1.0e-4)+1.0e-4
      q2=Ifi*IIKF+Iri*IIKR
      if(QBM.LT.0.5)then
         xvar2=1.0/NKF
         xvar3=q1**xvar2
         xvar1=(xvar3+4.0*q2)
         xvar4=xvar1**NKF
         qb=0.5*(q1+xvar4)
      else
         xvar1=(1.0+4.0*q2)
         xvar2=xvar1**NKF
         qb=0.5*q1*(1.0+xvar2)
      endif
      Itzr=Iri/qb
      Itzf=Ifi/qb
      if(ISP.GT.0.0)then
         argi=Vbep/(NFP*Vtv)
         expi=exp(argi)
         argx=Vbci/(NFP*Vtv)
         expx=exp(argx)
         Ifp=ISP*(WSP*expi+(1.0-WSP)*expx-1.0)
         q2p=Ifp*IIKP
         qbp=0.5*(1.0+sqrt(1.0+4.0*q2p))
         argi=Vbcp/(NFP*Vtv)
         expi=exp(argi)
         Irp=ISP*(expi-1.0)
         Iccp=(Ifp-Irp)/qbp
      else
         Ifp=0.0
         qbp=1.0
         Iccp=0.0
      endif
      if(WBE.EQ.1.0)then
         argi=Vbei/(NEI*Vtv)
         expi=exp(argi)
         argn=Vbei/(NEN*Vtv)
         expn=exp(argn)
         if(VBBE.GT.0.0)then
            argx=(-VBBE-Vbei)/(NBBE*Vtv)
            expx=exp(argx)
            Ibe=IBEI*(expi-1.0)+IBEN*(expn-1.0)-IBBE*(expx-EBBE)
         else
            Ibe=IBEI*(expi-1.0)+IBEN*(expn-1.0)
         endif
         Ibex=0.0
      elseif(WBE.EQ.0.0)then
         Ibe=0.0
         argi=Vbex/(NEI*Vtv)
         expi=exp(argi)
         argn=Vbex/(NEN*Vtv)
         expn=exp(argn)
         if(VBBE.GT.0.0)then
            argx=(-VBBE-Vbex)/(NBBE*Vtv)
            expx=exp(argx)
            Ibex=IBEI*(expi-1.0)+IBEN*(expn-1.0)-IBBE*(expx-EBBE)
         else
            Ibex=IBEI*(expi-1.0)+IBEN*(expn-1.0)
         endif
      else
         argi=Vbei/(NEI*Vtv)
         expi=exp(argi)
         argn=Vbei/(NEN*Vtv)
         expn=exp(argn)
         if(VBBE.GT.0.0)then
            argx=(-VBBE-Vbei)/(NBBE*Vtv)
            expx=exp(argx)
            Ibe=WBE*(IBEI*(expi-1.0)+IBEN*(expn-1.0)-IBBE*(expx-EBBE))
         else
            Ibe=WBE*(IBEI*(expi-1.0)+IBEN*(expn-1.0))
         endif
         argi=Vbex/(NEI*Vtv)
         expi=exp(argi)
         argn=Vbex/(NEN*Vtv)
         expn=exp(argn)
         if(VBBE.GT.0.0)then
            argx=(-VBBE-Vbex)/(NBBE*Vtv)
            expx=exp(argx)
            Ibex=(1.0-WBE)*(IBEI*(expi-1.0)+IBEN*(expn-1.0)-IBBE*(expx-E
     +BBE))
         else
            Ibex=(1.0-WBE)*(IBEI*(expi-1.0)+IBEN*(expn-1.0))
         endif
      endif
      argi=Vbci/(NCI*Vtv)
      expi=exp(argi)
      argn=Vbci/(NCN*Vtv)
      expn=exp(argn)
      Ibcj=IBCI*(expi-1.0)+IBCN*(expn-1.0)
      if((IBEIP.GT.0.0).OR.(IBENP.GT.0.0))then
         argi=Vbep/(NCI*Vtv)
         expi=exp(argi)
         argn=Vbep/(NCN*Vtv)
         expn=exp(argn)
         Ibep=IBEIP*(expi-1.0)+IBENP*(expn-1.0)
      else
         Ibep=0.0
      endif
      if(AVC1.GT.0.0)then
         vl=0.5*(sqrt((PC-Vbci)*(PC-Vbci)+0.01)+(PC-Vbci))
         xvar2=(MC-1.0)
         xvar3=vl**xvar2
         xvar1=-AVC2*xvar3
         xvar4=exp(xvar1)
         avalf=AVC1*vl*xvar4
         Igc=(Itzf-Itzr-Ibcj)*avalf
      else
         Igc=0.0
      endif
      Ibc=Ibcj-Igc
      if(RCX.GT.0.0)then
         Ircx=Vrcx/RCX
      else
         Ircx=0.0
      endif
      argi=Vbci/Vtv
      expi=exp(argi)
      argx=Vbcx/Vtv
      expx=exp(argx)
      Kbci=sqrt(1.0+GAMM*expi)
      Kbcx=sqrt(1.0+GAMM*expx)
      if(RCI.GT.0.0)then
         rKp1=(Kbci+1.0)/(Kbcx+1.0)
         xvar1=log(rKp1)
         Iohm=(Vrci+Vtv*(Kbci-Kbcx-xvar1))/RCI
         derf=IVO*RCI*Iohm/(1.0+0.5*IVO*IHRCF*sqrt(Vrci*Vrci+0.01))
         Irci=Iohm/sqrt(1.0+derf*derf)
      else
         Irci=0.0
      endif
      if(RBX.GT.0.0)then
         Irbx=Vrbx/RBX
      else
         Irbx=0.0
      endif
      if(RBI.GT.0.0)then
         Irbi=Vrbi*qb/RBI
      else
         Irbi=0.0
      endif
      if(RE.GT.0.0)then
         Ire=Vre/RE
      else
         Ire=0.0
      endif
      if(RBP.GT.0.0)then
         Irbp=Vrbp*qbp/RBP
      else
         Irbp=0.0
      endif
      if((IBCIP.GT.0.0).OR.(IBCNP.GT.0.0))then
         argi=Vbcp/(NCIP*Vtv)
         expi=exp(argi)
         argn=Vbcp/(NCNP*Vtv)
         expn=exp(argn)
         Ibcp=IBCIP*(expi-1.0)+IBCNP*(expn-1.0)
      else
         Ibcp=0.0
      endif
      if(RS.GT.0.0)then
         Irs=Vrs/RS
      else
         Irs=0.0
      endif
      if(Ifi.GT.0.0)then
         sgIf=1.0
      else
         sgIf=0.0
      endif
      rIf=Ifi*sgIf*IITF
      mIf=rIf/(rIf+1.0)
      xvar1=Vbci*IVTF/1.44
      xvar2=exp(xvar1)
      tff=TF*(1.0+QTF*q1)*(1.0+XTF*xvar2*(slTF+mIf*mIf)*sgIf)
      Qbe=CJE*qdbe*WBE+tff*Ifi/qb
      Qbex=CJE*qdbex*(1.0-WBE)
      Qbc=CJC*qdbc+TR*Iri+QCO*Kbci
      Qbcx=QCO*Kbcx
      Qbep=CJEP*qdbep+TR*Ifp
      Qbcp=CJCP*qdbcp+CCSO*Vbcp
      Qbeo=Vbe*CBEO
      Qbco=Vbc*CBCO
c
c     Scale outputs
c
      if(SCALE.NE.1.0)then
         Ibe=SCALE*Ibe
         Ibex=SCALE*Ibex
         Itzf=SCALE*Itzf
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
         Ibcp=SCALE*Ibcp
         Iccp=SCALE*Iccp
         Irs=SCALE*Irs
         Qbcp=SCALE*Qbcp
      endif
      return
      end
