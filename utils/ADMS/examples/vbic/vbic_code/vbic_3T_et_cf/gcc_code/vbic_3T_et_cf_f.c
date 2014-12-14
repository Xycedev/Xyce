
#include <stdio.h>
#include <math.h>

int vbic_3T_et_cf_f(p
	,Vrth,Vbei,Vbex,Vbci,Vbep,Vrcx
	,Vbcx,Vrci,Vrbx,Vrbi,Vre,Vrbp,Vbe
	,Vbc,Vcei,Ibe,Ibex,Itzf,Itzr,Ibc
	,Ibep,Ircx,Irci,Irbx,Irbi,Ire,Irbp
	,Qbe,Qbex,Qbc,Qbcx,Qbep,Qbeo,Qbco
	,Irth,Ith,Qcth,SCALE)
double *p
	,*Vrth,*Vbei,*Vbex,*Vbci,*Vbep,*Vrcx
	,*Vbcx,*Vrci,*Vrbx,*Vrbi,*Vre,*Vrbp,*Vbe
	,*Vbc,*Vcei,*Ibe,*Ibex,*Itzf,*Itzr,*Ibc
	,*Ibep,*Ircx,*Irci,*Irbx,*Irbi,*Ire,*Irbp
	,*Qbe,*Qbex,*Qbc,*Qbcx,*Qbep,*Qbeo,*Qbco
	,*Irth,*Ith,*Qcth,*SCALE;
{
double	Tini,Tdev,Vtv,rT,dT,xvar1,IKFatT;
double	RCXatT,RCIatT,RBXatT,RBIatT,REatT,RSatT,RBPatT;
double	xvar2,xvar3,xvar4,xvar5,xvar6,ISatT,ISRRatT;
double	ISPatT,IBEIatT,IBENatT,IBCIatT,IBCNatT,IBEIPatT,IBENPatT;
double	IBCIPatT,IBCNPatT,NFatT,NRatT,AVC2atT,VBBEatT,NBBEatT;
double	psiio,psiin,PEatT,PCatT,PSatT,CJEatT,CJCatT;
double	CJEPatT,CJCPatT,GAMMatT,VOatT,EBBEatT,IVEF,IVER;
double	IIKF,IIKR,IIKP,IVO,IHRCF,IVTF,IITF;
double	slTF,dv0,dvh,pwq,qlo,qhi,qdbe;
double	mv0,vl0,q0,dv,mv,vl,qdbex;
double	qdbc,vn0,vnl0,qlo0,vn,vnl,sel;
double	crt,cmx,cl,ql,qdbep,argi,expi;
double	Ifi,Iri,q1z,q1,q2,qb,argx;
double	expx,Ifp,q2p,qbp,argn,expn,Ibcj;
double	avalf,Igc,Kbci,Kbcx,rKp1,Iohm,derf;
double	sgIf,rIf,mIf,tff;

/*	Function only code */

	Tini=2.731500e+02+p[0];
	Tdev=(2.731500e+02+p[0])+(*Vrth);
	Vtv=1.380662e-23*Tdev/1.602189e-19;
	rT=Tdev/Tini;
	dT=Tdev-Tini;
	xvar1=pow(rT,p[90]);
	IKFatT=p[53]*xvar1;
	xvar1=pow(rT,p[91]);
	RCXatT=p[1]*xvar1;
	xvar1=pow(rT,p[68]);
	RCIatT=p[2]*xvar1;
	xvar1=pow(rT,p[92]);
	RBXatT=p[6]*xvar1;
	xvar1=pow(rT,p[67]);
	RBIatT=p[7]*xvar1;
	xvar1=pow(rT,p[66]);
	REatT=p[8]*xvar1;
	xvar1=pow(rT,p[69]);
	RSatT=p[9]*xvar1;
	xvar1=pow(rT,p[93]);
	RBPatT=p[10]*xvar1;
	xvar2=pow(rT,p[78]);
	xvar3=-p[71]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[12]);
	xvar6=pow(xvar1,xvar5);
	ISatT=p[11]*xvar6;
	xvar2=pow(rT,p[95]);
	xvar3=-p[96]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[13]);
	xvar6=pow(xvar1,xvar5);
	ISRRatT=p[94]*xvar6;
	xvar2=pow(rT,p[78]);
	xvar3=-p[97]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[44]);
	xvar6=pow(xvar1,xvar5);
	ISPatT=p[42]*xvar6;
	xvar2=pow(rT,p[79]);
	xvar3=-p[72]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[33]);
	xvar6=pow(xvar1,xvar5);
	IBEIatT=p[31]*xvar6;
	xvar2=pow(rT,p[80]);
	xvar3=-p[75]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[35]);
	xvar6=pow(xvar1,xvar5);
	IBENatT=p[34]*xvar6;
	xvar2=pow(rT,p[79]);
	xvar3=-p[73]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[37]);
	xvar6=pow(xvar1,xvar5);
	IBCIatT=p[36]*xvar6;
	xvar2=pow(rT,p[80]);
	xvar3=-p[76]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[39]);
	xvar6=pow(xvar1,xvar5);
	IBCNatT=p[38]*xvar6;
	xvar2=pow(rT,p[79]);
	xvar3=-p[73]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[37]);
	xvar6=pow(xvar1,xvar5);
	IBEIPatT=p[45]*xvar6;
	xvar2=pow(rT,p[80]);
	xvar3=-p[76]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[39]);
	xvar6=pow(xvar1,xvar5);
	IBENPatT=p[46]*xvar6;
	xvar2=pow(rT,p[79]);
	xvar3=-p[74]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[48]);
	xvar6=pow(xvar1,xvar5);
	IBCIPatT=p[47]*xvar6;
	xvar2=pow(rT,p[80]);
	xvar3=-p[77]*(1.0-rT)/Vtv;
	xvar4=exp(xvar3);
	xvar1=(xvar2*xvar4);
	xvar5=(1.0/p[50]);
	xvar6=pow(xvar1,xvar5);
	IBCNPatT=p[49]*xvar6;
	NFatT=p[12]*(1.0+dT*p[81]);
	NRatT=p[13]*(1.0+dT*p[81]);
	AVC2atT=p[41]*(1.0+dT*p[82]);
	VBBEatT=p[98]*(1.0+dT*(p[101]+dT*p[102]));
	NBBEatT=p[99]*(1.0+dT*p[103]);
	xvar2=0.5*p[17]*rT/Vtv;
	xvar3=exp(xvar2);
	xvar4=-0.5*p[17]*rT/Vtv;
	xvar5=exp(xvar4);
	xvar1=xvar3-xvar5;
	xvar6=log(xvar1);
	psiio=2.0*(Vtv/rT)*xvar6;
	xvar1=log(rT);
	psiin=psiio*rT-3.0*Vtv*xvar1-p[72]*(rT-1.0);
	xvar2=-psiin/Vtv;
	xvar3=exp(xvar2);
	xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3));
	xvar4=log(xvar1);
	PEatT=psiin+2.0*Vtv*xvar4;
	xvar2=0.5*p[24]*rT/Vtv;
	xvar3=exp(xvar2);
	xvar4=-0.5*p[24]*rT/Vtv;
	xvar5=exp(xvar4);
	xvar1=xvar3-xvar5;
	xvar6=log(xvar1);
	psiio=2.0*(Vtv/rT)*xvar6;
	xvar1=log(rT);
	psiin=psiio*rT-3.0*Vtv*xvar1-p[73]*(rT-1.0);
	xvar2=-psiin/Vtv;
	xvar3=exp(xvar2);
	xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3));
	xvar4=log(xvar1);
	PCatT=psiin+2.0*Vtv*xvar4;
	xvar2=0.5*p[28]*rT/Vtv;
	xvar3=exp(xvar2);
	xvar4=-0.5*p[28]*rT/Vtv;
	xvar5=exp(xvar4);
	xvar1=xvar3-xvar5;
	xvar6=log(xvar1);
	psiio=2.0*(Vtv/rT)*xvar6;
	xvar1=log(rT);
	psiin=psiio*rT-3.0*Vtv*xvar1-p[74]*(rT-1.0);
	xvar2=-psiin/Vtv;
	xvar3=exp(xvar2);
	xvar1=0.5*(1.0+sqrt(1.0+4.0*xvar3));
	xvar4=log(xvar1);
	PSatT=psiin+2.0*Vtv*xvar4;
	xvar1=p[17]/PEatT;
	xvar2=pow(xvar1,p[18]);
	CJEatT=p[16]*xvar2;
	xvar1=p[24]/PCatT;
	xvar2=pow(xvar1,p[25]);
	CJCatT=p[21]*xvar2;
	xvar1=p[24]/PCatT;
	xvar2=pow(xvar1,p[25]);
	CJEPatT=p[23]*xvar2;
	xvar1=p[28]/PSatT;
	xvar2=pow(xvar1,p[29]);
	CJCPatT=p[27]*xvar2;
	xvar1=pow(rT,p[78]);
	xvar2=-p[71]*(1.0-rT)/Vtv;
	xvar3=exp(xvar2);
	GAMMatT=p[4]*xvar1*xvar3;
	xvar1=pow(rT,p[70]);
	VOatT=p[3]*xvar1;
	xvar1=-VBBEatT/(NBBEatT*Vtv);
	EBBEatT=exp(xvar1);
	if(p[51]>0.0){
		IVEF=1.0/p[51];
	}else{
		IVEF=0.0;
	}
	if(p[52]>0.0){
		IVER=1.0/p[52];
	}else{
		IVER=0.0;
	}
	if(p[53]>0.0){
		IIKF=1.0/IKFatT;
	}else{
		IIKF=0.0;
	}
	if(p[54]>0.0){
		IIKR=1.0/p[54];
	}else{
		IIKR=0.0;
	}
	if(p[55]>0.0){
		IIKP=1.0/p[55];
	}else{
		IIKP=0.0;
	}
	if(p[3]>0.0){
		IVO=1.0/VOatT;
	}else{
		IVO=0.0;
	}
	if(p[5]>0.0){
		IHRCF=1.0/p[5];
	}else{
		IHRCF=0.0;
	}
	if(p[59]>0.0){
		IVTF=1.0/p[59];
	}else{
		IVTF=0.0;
	}
	if(p[60]>0.0){
		IITF=1.0/p[60];
	}else{
		IITF=0.0;
	}
	if(p[60]>0.0){
		slTF=0.0;
	}else{
		slTF=1.0;
	}
	dv0=-PEatT*p[14];
	if(p[19]<=0.0){
		dvh=(*Vbei)+dv0;
		if(dvh>0.0){
			xvar1=(1.0-p[14]);
			xvar2=(-1.0-p[18]);
			pwq=pow(xvar1,xvar2);
			qlo=PEatT*(1.0-pwq*(1.0-p[14])*(1.0-p[14]))/(1.0-p[18]);
			qhi=dvh*(1.0-p[14]+0.5*p[18]*dvh/PEatT)*pwq;
		}else{
			xvar1=(1.0-(*Vbei)/PEatT);
			xvar2=(1.0-p[18]);
			xvar3=pow(xvar1,xvar2);
			qlo=PEatT*(1.0-xvar3)/(1.0-p[18]);
			qhi=0.0;
		}
		qdbe=qlo+qhi;
	}else{
		mv0=sqrt(dv0*dv0+4.0*p[19]*p[19]);
		vl0=-0.5*(dv0+mv0);
		xvar1=(1.0-vl0/PEatT);
		xvar2=(1.0-p[18]);
		xvar3=pow(xvar1,xvar2);
		q0=-PEatT*xvar3/(1.0-p[18]);
		dv=(*Vbei)+dv0;
		mv=sqrt(dv*dv+4.0*p[19]*p[19]);
		vl=0.5*(dv-mv)-dv0;
		xvar1=(1.0-vl/PEatT);
		xvar2=(1.0-p[18]);
		xvar3=pow(xvar1,xvar2);
		qlo=-PEatT*xvar3/(1.0-p[18]);
		xvar1=(1.0-p[14]);
		xvar2=(-p[18]);
		xvar3=pow(xvar1,xvar2);
		qdbe=qlo+xvar3*((*Vbei)-vl+vl0)-q0;
	}
	dv0=-PEatT*p[14];
	if(p[19]<=0.0){
		dvh=(*Vbex)+dv0;
		if(dvh>0.0){
			xvar1=(1.0-p[14]);
			xvar2=(-1.0-p[18]);
			pwq=pow(xvar1,xvar2);
			qlo=PEatT*(1.0-pwq*(1.0-p[14])*(1.0-p[14]))/(1.0-p[18]);
			qhi=dvh*(1.0-p[14]+0.5*p[18]*dvh/PEatT)*pwq;
		}else{
			xvar1=(1.0-(*Vbex)/PEatT);
			xvar2=(1.0-p[18]);
			xvar3=pow(xvar1,xvar2);
			qlo=PEatT*(1.0-xvar3)/(1.0-p[18]);
			qhi=0.0;
		}
		qdbex=qlo+qhi;
	}else{
		mv0=sqrt(dv0*dv0+4.0*p[19]*p[19]);
		vl0=-0.5*(dv0+mv0);
		xvar1=(1.0-vl0/PEatT);
		xvar2=(1.0-p[18]);
		xvar3=pow(xvar1,xvar2);
		q0=-PEatT*xvar3/(1.0-p[18]);
		dv=(*Vbex)+dv0;
		mv=sqrt(dv*dv+4.0*p[19]*p[19]);
		vl=0.5*(dv-mv)-dv0;
		xvar1=(1.0-vl/PEatT);
		xvar2=(1.0-p[18]);
		xvar3=pow(xvar1,xvar2);
		qlo=-PEatT*xvar3/(1.0-p[18]);
		xvar1=(1.0-p[14]);
		xvar2=(-p[18]);
		xvar3=pow(xvar1,xvar2);
		qdbex=qlo+xvar3*((*Vbex)-vl+vl0)-q0;
	}
	dv0=-PCatT*p[14];
	if(p[26]<=0.0){
		dvh=(*Vbci)+dv0;
		if(dvh>0.0){
			xvar1=(1.0-p[14]);
			xvar2=(-1.0-p[25]);
			pwq=pow(xvar1,xvar2);
			qlo=PCatT*(1.0-pwq*(1.0-p[14])*(1.0-p[14]))/(1.0-p[25]);
			qhi=dvh*(1.0-p[14]+0.5*p[25]*dvh/PCatT)*pwq;
		}else{
			if((p[85]>0.0)&&((*Vbci)<-p[85])){
				xvar1=(1.0+p[85]/PCatT);
				xvar2=(1.0-p[25]);
				xvar3=pow(xvar1,xvar2);
				qlo=PCatT*(1.0-xvar3*(1.0-((1.0-p[25])*((*Vbci)+p[85]))/(PCatT+p[85])))/(1.0-p[25]);
			}else{
				xvar1=(1.0-(*Vbci)/PCatT);
				xvar2=(1.0-p[25]);
				xvar3=pow(xvar1,xvar2);
				qlo=PCatT*(1.0-xvar3)/(1.0-p[25]);
			}
			qhi=0.0;
		}
		qdbc=qlo+qhi;
	}else{
		if((p[85]>0.0)&&(p[86]>0.0)){
			vn0=(p[85]+dv0)/(p[85]-dv0);
			vnl0=2.0*vn0/(sqrt((vn0-1.0)*(vn0-1.0)+4.0*p[26]*p[26])+sqrt((vn0+1.0)*(vn0+1.0)+4.0*p[86]*p[86]));
			vl0=0.5*(vnl0*(p[85]-dv0)-p[85]-dv0);
			xvar1=(1.0-vl0/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			qlo0=PCatT*(1.0-xvar3)/(1.0-p[25]);
			vn=(2.0*(*Vbci)+p[85]+dv0)/(p[85]-dv0);
			vnl=2.0*vn/(sqrt((vn-1.0)*(vn-1.0)+4.0*p[26]*p[26])+sqrt((vn+1.0)*(vn+1.0)+4.0*p[86]*p[86]));
			vl=0.5*(vnl*(p[85]-dv0)-p[85]-dv0);
			xvar1=(1.0-vl/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			qlo=PCatT*(1.0-xvar3)/(1.0-p[25]);
			sel=0.5*(vnl+1.0);
			xvar1=(1.0+p[85]/PCatT);
			xvar2=(-p[25]);
			crt=pow(xvar1,xvar2);
			xvar1=(1.0+dv0/PCatT);
			xvar2=(-p[25]);
			cmx=pow(xvar1,xvar2);
			cl=(1.0-sel)*crt+sel*cmx;
			ql=((*Vbci)-vl+vl0)*cl;
			qdbc=ql+qlo-qlo0;
		}else{
			mv0=sqrt(dv0*dv0+4.0*p[26]*p[26]);
			vl0=-0.5*(dv0+mv0);
			xvar1=(1.0-vl0/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			q0=-PCatT*xvar3/(1.0-p[25]);
			dv=(*Vbci)+dv0;
			mv=sqrt(dv*dv+4.0*p[26]*p[26]);
			vl=0.5*(dv-mv)-dv0;
			xvar1=(1.0-vl/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			qlo=-PCatT*xvar3/(1.0-p[25]);
			xvar1=(1.0-p[14]);
			xvar2=(-p[25]);
			xvar3=pow(xvar1,xvar2);
			qdbc=qlo+xvar3*((*Vbci)-vl+vl0)-q0;
		}
	}
	dv0=-PCatT*p[14];
	if(p[26]<=0.0){
		dvh=(*Vbep)+dv0;
		if(dvh>0.0){
			xvar1=(1.0-p[14]);
			xvar2=(-1.0-p[25]);
			pwq=pow(xvar1,xvar2);
			qlo=PCatT*(1.0-pwq*(1.0-p[14])*(1.0-p[14]))/(1.0-p[25]);
			qhi=dvh*(1.0-p[14]+0.5*p[25]*dvh/PCatT)*pwq;
		}else{
			if((p[85]>0.0)&&((*Vbep)<-p[85])){
				xvar1=(1.0+p[85]/PCatT);
				xvar2=(1.0-p[25]);
				xvar3=pow(xvar1,xvar2);
				qlo=PCatT*(1.0-xvar3*(1.0-((1.0-p[25])*((*Vbep)+p[85]))/(PCatT+p[85])))/(1.0-p[25]);
			}else{
				xvar1=(1.0-(*Vbep)/PCatT);
				xvar2=(1.0-p[25]);
				xvar3=pow(xvar1,xvar2);
				qlo=PCatT*(1.0-xvar3)/(1.0-p[25]);
			}
			qhi=0.0;
		}
		qdbep=qlo+qhi;
	}else{
		if((p[85]>0.0)&&(p[86]>0.0)){
			vn0=(p[85]+dv0)/(p[85]-dv0);
			vnl0=2.0*vn0/(sqrt((vn0-1.0)*(vn0-1.0)+4.0*p[26]*p[26])+sqrt((vn0+1.0)*(vn0+1.0)+4.0*p[86]*p[86]));
			vl0=0.5*(vnl0*(p[85]-dv0)-p[85]-dv0);
			xvar1=(1.0-vl0/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			qlo0=PCatT*(1.0-xvar3)/(1.0-p[25]);
			vn=(2.0*(*Vbep)+p[85]+dv0)/(p[85]-dv0);
			vnl=2.0*vn/(sqrt((vn-1.0)*(vn-1.0)+4.0*p[26]*p[26])+sqrt((vn+1.0)*(vn+1.0)+4.0*p[86]*p[86]));
			vl=0.5*(vnl*(p[85]-dv0)-p[85]-dv0);
			xvar1=(1.0-vl/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			qlo=PCatT*(1.0-xvar3)/(1.0-p[25]);
			sel=0.5*(vnl+1.0);
			xvar1=(1.0+p[85]/PCatT);
			xvar2=(-p[25]);
			crt=pow(xvar1,xvar2);
			xvar1=(1.0+dv0/PCatT);
			xvar2=(-p[25]);
			cmx=pow(xvar1,xvar2);
			cl=(1.0-sel)*crt+sel*cmx;
			ql=((*Vbep)-vl+vl0)*cl;
			qdbep=ql+qlo-qlo0;
		}else{
			mv0=sqrt(dv0*dv0+4.0*p[26]*p[26]);
			vl0=-0.5*(dv0+mv0);
			xvar1=(1.0-vl0/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			q0=-PCatT*xvar3/(1.0-p[25]);
			dv=(*Vbep)+dv0;
			mv=sqrt(dv*dv+4.0*p[26]*p[26]);
			vl=0.5*(dv-mv)-dv0;
			xvar1=(1.0-vl/PCatT);
			xvar2=(1.0-p[25]);
			xvar3=pow(xvar1,xvar2);
			qlo=-PCatT*xvar3/(1.0-p[25]);
			xvar1=(1.0-p[14]);
			xvar2=(-p[25]);
			xvar3=pow(xvar1,xvar2);
			qdbep=qlo+xvar3*((*Vbep)-vl+vl0)-q0;
		}
	}
	argi=(*Vbei)/(NFatT*Vtv);
	expi=exp(argi);
	Ifi=ISatT*(expi-1.0);
	argi=(*Vbci)/(NRatT*Vtv);
	expi=exp(argi);
	Iri=ISatT*ISRRatT*(expi-1.0);
	q1z=1.0+qdbe*IVER+qdbc*IVEF;
	q1=0.5*(sqrt((q1z-1.0e-4)*(q1z-1.0e-4)+1.0e-8)+q1z-1.0e-4)+1.0e-4;
	q2=Ifi*IIKF+Iri*IIKR;
	if(p[88]<0.5){
		xvar2=1.0/p[89];
		xvar3=pow(q1,xvar2);
		xvar1=(xvar3+4.0*q2);
		xvar4=pow(xvar1,p[89]);
		qb=0.5*(q1+xvar4);
	}else{
		xvar1=(1.0+4.0*q2);
		xvar2=pow(xvar1,p[89]);
		qb=0.5*q1*(1.0+xvar2);
	}
	(*Itzr)=Iri/qb;
	(*Itzf)=Ifi/qb;
	if(p[42]>0.0){
		argi=(*Vbep)/(p[44]*Vtv);
		expi=exp(argi);
		argx=(*Vbci)/(p[44]*Vtv);
		expx=exp(argx);
		Ifp=ISPatT*(p[43]*expi+(1.0-p[43])*expx-1.0);
		q2p=Ifp*IIKP;
		qbp=0.5*(1.0+sqrt(1.0+4.0*q2p));
	}else{
		Ifp=0.0;
		qbp=1.0;
	}
	if(p[32]==1.0){
		argi=(*Vbei)/(p[33]*Vtv);
		expi=exp(argi);
		argn=(*Vbei)/(p[35]*Vtv);
		expn=exp(argn);
		if(p[98]>0.0){
			argx=(-VBBEatT-(*Vbei))/(NBBEatT*Vtv);
			expx=exp(argx);
			(*Ibe)=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p[100]*(expx-EBBEatT);
		}else{
			(*Ibe)=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0);
		}
		(*Ibex)=0.0;
	}else if(p[32]==0.0){
		(*Ibe)=0.0;
		argi=(*Vbex)/(p[33]*Vtv);
		expi=exp(argi);
		argn=(*Vbex)/(p[35]*Vtv);
		expn=exp(argn);
		if(p[98]>0.0){
			argx=(-VBBEatT-(*Vbex))/(NBBEatT*Vtv);
			expx=exp(argx);
			(*Ibex)=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p[100]*(expx-EBBEatT);
		}else{
			(*Ibex)=IBEIatT*(expi-1.0)+IBENatT*(expn-1.0);
		}
	}else{
		argi=(*Vbei)/(p[33]*Vtv);
		expi=exp(argi);
		argn=(*Vbei)/(p[35]*Vtv);
		expn=exp(argn);
		if(p[98]>0.0){
			argx=(-VBBEatT-(*Vbei))/(NBBEatT*Vtv);
			expx=exp(argx);
			(*Ibe)=p[32]*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p[100]*(expx-EBBEatT));
		}else{
			(*Ibe)=p[32]*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0));
		}
		argi=(*Vbex)/(p[33]*Vtv);
		expi=exp(argi);
		argn=(*Vbex)/(p[35]*Vtv);
		expn=exp(argn);
		if(p[98]>0.0){
			argx=(-VBBEatT-(*Vbex))/(NBBEatT*Vtv);
			expx=exp(argx);
			(*Ibex)=(1.0-p[32])*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0)-p[100]*(expx-EBBEatT));
		}else{
			(*Ibex)=(1.0-p[32])*(IBEIatT*(expi-1.0)+IBENatT*(expn-1.0));
		}
	}
	argi=(*Vbci)/(p[37]*Vtv);
	expi=exp(argi);
	argn=(*Vbci)/(p[39]*Vtv);
	expn=exp(argn);
	Ibcj=IBCIatT*(expi-1.0)+IBCNatT*(expn-1.0);
	if((p[45]>0.0)||(p[46]>0.0)){
		argi=(*Vbep)/(p[37]*Vtv);
		expi=exp(argi);
		argn=(*Vbep)/(p[39]*Vtv);
		expn=exp(argn);
		(*Ibep)=IBEIPatT*(expi-1.0)+IBENPatT*(expn-1.0);
	}else{
		(*Ibep)=0.0;
	}
	if(p[40]>0.0){
		vl=0.5*(sqrt((PCatT-(*Vbci))*(PCatT-(*Vbci))+0.01)+(PCatT-(*Vbci)));
		xvar2=(p[25]-1.0);
		xvar3=pow(vl,xvar2);
		xvar1=-AVC2atT*xvar3;
		xvar4=exp(xvar1);
		avalf=p[40]*vl*xvar4;
		Igc=((*Itzf)-(*Itzr)-Ibcj)*avalf;
	}else{
		Igc=0.0;
	}
	(*Ibc)=Ibcj-Igc;
	if(p[1]>0.0){
		(*Ircx)=(*Vrcx)/RCXatT;
	}else{
		(*Ircx)=0.0;
	}
	argi=(*Vbci)/Vtv;
	expi=exp(argi);
	argx=(*Vbcx)/Vtv;
	expx=exp(argx);
	Kbci=sqrt(1.0+GAMMatT*expi);
	Kbcx=sqrt(1.0+GAMMatT*expx);
	if(p[2]>0.0){
		rKp1=(Kbci+1.0)/(Kbcx+1.0);
		xvar1=log(rKp1);
		Iohm=((*Vrci)+Vtv*(Kbci-Kbcx-xvar1))/RCIatT;
		derf=IVO*RCIatT*Iohm/(1.0+0.5*IVO*IHRCF*sqrt((*Vrci)*(*Vrci)+0.01));
		(*Irci)=Iohm/sqrt(1.0+derf*derf);
	}else{
		(*Irci)=0.0;
	}
	if(p[6]>0.0){
		(*Irbx)=(*Vrbx)/RBXatT;
	}else{
		(*Irbx)=0.0;
	}
	if(p[7]>0.0){
		(*Irbi)=(*Vrbi)*qb/RBIatT;
	}else{
		(*Irbi)=0.0;
	}
	if(p[8]>0.0){
		(*Ire)=(*Vre)/REatT;
	}else{
		(*Ire)=0.0;
	}
	if(p[10]>0.0){
		(*Irbp)=(*Vrbp)*qbp/RBPatT;
	}else{
		(*Irbp)=0.0;
	}
	if(Ifi>0.0){
		sgIf=1.0;
	}else{
		sgIf=0.0;
	}
	rIf=Ifi*sgIf*IITF;
	mIf=rIf/(rIf+1.0);
	xvar1=(*Vbci)*IVTF/1.44;
	xvar2=exp(xvar1);
	tff=p[56]*(1.0+p[57]*q1)*(1.0+p[58]*xvar2*(slTF+mIf*mIf)*sgIf);
	(*Qbe)=CJEatT*qdbe*p[32]+tff*Ifi/qb;
	(*Qbex)=CJEatT*qdbex*(1.0-p[32]);
	(*Qbc)=CJCatT*qdbc+p[61]*Iri+p[22]*Kbci;
	(*Qbcx)=p[22]*Kbcx;
	(*Qbep)=CJEPatT*qdbep+p[61]*Ifp;
	(*Qbeo)=(*Vbe)*p[15];
	(*Qbco)=(*Vbc)*p[20];
	(*Ith)=-((*Ibe)*(*Vbei)+(*Ibc)*(*Vbci)+((*Itzf)-(*Itzr))*(*Vcei)+(*Ibex)*(*Vbex)+(*Ibep)*(*Vbep)+(*Ircx)*(*Vrcx)+(*Irci)*(*Vrci)+(*Irbx)*(*Vrbx)+(*Irbi)*(*Vrbi)+(*Ire)*(*Vre)+(*Irbp)*(*Vrbp));
	if(p[83]>0.0){
		(*Irth)=(*Vrth)/p[83];
	}else{
		(*Irth)=0.0;
	}
	(*Qcth)=(*Vrth)*p[84];

/*	Scale outputs */
	if(*SCALE!=1.0){
		*Ibe=(*SCALE)*(*Ibe);
		*Ibex=(*SCALE)*(*Ibex);
		*Itzf=(*SCALE)*(*Itzf);
		*Itzr=(*SCALE)*(*Itzr);
		*Ibc=(*SCALE)*(*Ibc);
		*Ibep=(*SCALE)*(*Ibep);
		*Ircx=(*SCALE)*(*Ircx);
		*Irci=(*SCALE)*(*Irci);
		*Irbx=(*SCALE)*(*Irbx);
		*Irbi=(*SCALE)*(*Irbi);
		*Ire=(*SCALE)*(*Ire);
		*Irbp=(*SCALE)*(*Irbp);
		*Qbe=(*SCALE)*(*Qbe);
		*Qbex=(*SCALE)*(*Qbex);
		*Qbc=(*SCALE)*(*Qbc);
		*Qbcx=(*SCALE)*(*Qbcx);
		*Qbep=(*SCALE)*(*Qbep);
		*Qbeo=(*SCALE)*(*Qbeo);
		*Qbco=(*SCALE)*(*Qbco);
		*Irth=(*SCALE)*(*Irth);
		*Ith=(*SCALE)*(*Ith);
		*Qcth=(*SCALE)*(*Qcth);
	}
	return(0);
}
