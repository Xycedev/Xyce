c
c	vbic_3T_et_xf Solver
c
c	This program is a stand-alone FORTRAN solver for vbic_3T_et_xf.
c	For usage see the help message that is printed below.
c
	implicit double precision (a-z)
	double precision p(108),pnom(108)
	double precision xx(3),cx(3),vmx(4),imx(3)
	double precision x(10),xlast(10),xinit(10)
	double precision G(11,11),C(11,11),xxlast(3)
	double precision Ginv(11,11),Glu(11,11),B(11,11)
	double precision BGinv(11,11),Geff(11,11)
	double precision Frb(11),Fib(11),Frc(11),Fic(11)
	double precision fr(11),fi(11)
	integer i,j,k,nit,nac,iargc,nargs,acc,ierr,iinit,doac,dodc
	integer ipvt(11),ipvtGlu(11)
	character*256 arg,fname,pname
c
c	Bias specification node order is: c b e
c
	acc=0
	doac=0
	nac=11
c
c	Default parameter values
c
	pnom(1)=27.0
	pnom(2)=0.0
	pnom(3)=0.0
	pnom(4)=0.0
	pnom(5)=0.0
	pnom(6)=0.0
	pnom(7)=0.0
	pnom(8)=0.0
	pnom(9)=0.0
	pnom(10)=0.0
	pnom(11)=0.0
	pnom(12)=1.0e-16
	pnom(13)=1.0
	pnom(14)=1.0
	pnom(15)=0.9
	pnom(16)=0.0
	pnom(17)=0.0
	pnom(18)=0.75
	pnom(19)=0.33
	pnom(20)=-0.5
	pnom(21)=0.0
	pnom(22)=0.0
	pnom(23)=0.0
	pnom(24)=0.0
	pnom(25)=0.75
	pnom(26)=0.33
	pnom(27)=-0.5
	pnom(28)=0.0
	pnom(29)=0.75
	pnom(30)=0.33
	pnom(31)=-0.5
	pnom(32)=1.0e-18
	pnom(33)=1.0
	pnom(34)=1.0
	pnom(35)=0.0
	pnom(36)=2.0
	pnom(37)=1.0e-16
	pnom(38)=1.0
	pnom(39)=0.0
	pnom(40)=2.0
	pnom(41)=0.0
	pnom(42)=0.0
	pnom(43)=0.0
	pnom(44)=1.0
	pnom(45)=1.0
	pnom(46)=0.0
	pnom(47)=0.0
	pnom(48)=0.0
	pnom(49)=1.0
	pnom(50)=0.0
	pnom(51)=2.0
	pnom(52)=0.0
	pnom(53)=0.0
	pnom(54)=0.0
	pnom(55)=0.0
	pnom(56)=0.0
	pnom(57)=0.0
	pnom(58)=0.0
	pnom(59)=0.0
	pnom(60)=0.0
	pnom(61)=0.0
	pnom(62)=0.0
	pnom(63)=0.0
	pnom(64)=0.0
	pnom(65)=1.0
	pnom(66)=1.0
	pnom(67)=0
	pnom(68)=0
	pnom(69)=0
	pnom(70)=0
	pnom(71)=0
	pnom(72)=1.12
	pnom(73)=1.12
	pnom(74)=1.12
	pnom(75)=1.12
	pnom(76)=1.12
	pnom(77)=1.12
	pnom(78)=1.12
	pnom(79)=3.0
	pnom(80)=3.0
	pnom(81)=3.0
	pnom(82)=0.0
	pnom(83)=0.0
	pnom(84)=0.0
	pnom(85)=0.0
	pnom(86)=0.0
	pnom(87)=0.1
	pnom(88)=0.0
	pnom(89)=0.0
	pnom(90)=0.5
	pnom(91)=0
	pnom(92)=0
	pnom(93)=0
	pnom(94)=0
	pnom(95)=1.0
	pnom(96)=0.0
	pnom(97)=0.0
	pnom(98)=1.12
	pnom(99)=0.0
	pnom(100)=1.0
	pnom(101)=1.0e-6
	pnom(102)=0.0
	pnom(103)=0.0
	pnom(104)=0.0
	pnom(105)=0.0
	pnom(106)=0.0
	pnom(107)=1.2
	pnom(108)=0.0
c
c	Parse command line arguments
c
	TAMB=27.0d0
	fname='/none'
	nargs=iargc()
	i=1
20	continue
	if(i.GT.nargs)goto50
	call getarg(i,arg)
	i=i+1
	if(arg.EQ.'-a')then
		acc=1
		goto20
	endif
	if(arg.EQ.'-ac')then
		doac=1
		goto20
	endif
	if(arg.EQ.'-yp')then
		doac=1
		goto20
	endif
	if(arg.EQ.'-h')then
		write(6,*) ''
		write(6,*) 'vbic_3T_et_xf_solveri: solver for model vbic_3T_et_xf'
		write(6,*) ''
		write(6,*) 'Usage: vbic_3T_et_xf_solveri [options] [ParameterFileName]'
		write(6,*) '       followed by lines specifying terminal biases'
		write(6,*) ''
		write(6,*) 'Files:'
		write(6,*) '	ParameterFileName	File with parameter name/value pairs'
		write(6,*) '				(format is one name and value pair'
		write(6,*) '				 per line, default values are the'
		write(6,*) '				 vbic_3T_et_xf defaults)'
		write(6,*) ''
		write(6,*) 'Options:'
		write(6,*) '	-a			Write solution in 30.20e format'
		write(6,*) '	-h			Print this help message'
		write(6,*) '	-t TEMPC		Simulate at TEMP C (default is TNOM)'
		write(6,*) '	-yp			Do y-parameters'
		write(6,*) ''
		write(6,*) 'The terminal biases must all be specified on a single line.'
		write(6,*) 'Units are volts and amps and order is:'
		write(6,*) '	c b e dt [f for -yp option]'
		write(6,*) 'For vbic_3T_et_xf_solveri bias is current for the base node'
		write(6,*) 'The vbic_3T_et_xf solution is then printed, the first row of'
		write(6,*) 'results is the terminal voltages and the second row'
		write(6,*) 'is the terminal currents.'
		write(6,*) 'When specifying applied terminal biases from the keyboard'
		write(6,*) 'rather than via STDIN redirection from a file, enter Ctrl-D'
		write(6,*) '(EOF) to gracefully terminate this program.'
		write(6,*) ''
		call exit(0)
	endif
	if(arg.EQ.'-t')then
		if(i.GT.nargs)then
			write(6,*) 'ERROR: no temperature specified for -t flag'
			call exit(1)
		endif
		call getarg(i,arg)
		i=i+1
		read(arg,*,err=49) TAMB
		if(TAMB<=-273)goto49
		goto20
49		write(6,*) 'ERROR: TAMB must be a valid temperature'
		call exit(1)
	endif
	if(arg(1:1).EQ.'-')then
		write(6,*) 'ERROR: unknown flag: ',arg
		call exit(1)
	endif
	fname=arg
	goto20
50	continue
c
c	Read parameters from file
c
	if(fname.NE.'/none')then
		open(31,file=fname,err=57)
51		read(31,*,end=58) pname,pvalue
		if(pname.EQ.'TNOM')pnom(1)=pvalue
		if(pname.EQ.'RCX')pnom(2)=pvalue
		if(pname.EQ.'RCI')pnom(3)=pvalue
		if(pname.EQ.'VO')pnom(4)=pvalue
		if(pname.EQ.'GAMM')pnom(5)=pvalue
		if(pname.EQ.'HRCF')pnom(6)=pvalue
		if(pname.EQ.'RBX')pnom(7)=pvalue
		if(pname.EQ.'RBI')pnom(8)=pvalue
		if(pname.EQ.'RE')pnom(9)=pvalue
		if(pname.EQ.'RS')pnom(10)=pvalue
		if(pname.EQ.'RBP')pnom(11)=pvalue
		if(pname.EQ.'IS')pnom(12)=pvalue
		if(pname.EQ.'NF')pnom(13)=pvalue
		if(pname.EQ.'NR')pnom(14)=pvalue
		if(pname.EQ.'FC')pnom(15)=pvalue
		if(pname.EQ.'CBEO')pnom(16)=pvalue
		if(pname.EQ.'CJE')pnom(17)=pvalue
		if(pname.EQ.'PE')pnom(18)=pvalue
		if(pname.EQ.'ME')pnom(19)=pvalue
		if(pname.EQ.'AJE')pnom(20)=pvalue
		if(pname.EQ.'CBCO')pnom(21)=pvalue
		if(pname.EQ.'CJC')pnom(22)=pvalue
		if(pname.EQ.'QCO')pnom(23)=pvalue
		if(pname.EQ.'CJEP')pnom(24)=pvalue
		if(pname.EQ.'PC')pnom(25)=pvalue
		if(pname.EQ.'MC')pnom(26)=pvalue
		if(pname.EQ.'AJC')pnom(27)=pvalue
		if(pname.EQ.'CJCP')pnom(28)=pvalue
		if(pname.EQ.'PS')pnom(29)=pvalue
		if(pname.EQ.'MS')pnom(30)=pvalue
		if(pname.EQ.'AJS')pnom(31)=pvalue
		if(pname.EQ.'IBEI')pnom(32)=pvalue
		if(pname.EQ.'WBE')pnom(33)=pvalue
		if(pname.EQ.'NEI')pnom(34)=pvalue
		if(pname.EQ.'IBEN')pnom(35)=pvalue
		if(pname.EQ.'NEN')pnom(36)=pvalue
		if(pname.EQ.'IBCI')pnom(37)=pvalue
		if(pname.EQ.'NCI')pnom(38)=pvalue
		if(pname.EQ.'IBCN')pnom(39)=pvalue
		if(pname.EQ.'NCN')pnom(40)=pvalue
		if(pname.EQ.'AVC1')pnom(41)=pvalue
		if(pname.EQ.'AVC2')pnom(42)=pvalue
		if(pname.EQ.'ISP')pnom(43)=pvalue
		if(pname.EQ.'WSP')pnom(44)=pvalue
		if(pname.EQ.'NFP')pnom(45)=pvalue
		if(pname.EQ.'IBEIP')pnom(46)=pvalue
		if(pname.EQ.'IBENP')pnom(47)=pvalue
		if(pname.EQ.'IBCIP')pnom(48)=pvalue
		if(pname.EQ.'NCIP')pnom(49)=pvalue
		if(pname.EQ.'IBCNP')pnom(50)=pvalue
		if(pname.EQ.'NCNP')pnom(51)=pvalue
		if(pname.EQ.'VEF')pnom(52)=pvalue
		if(pname.EQ.'VER')pnom(53)=pvalue
		if(pname.EQ.'IKF')pnom(54)=pvalue
		if(pname.EQ.'IKR')pnom(55)=pvalue
		if(pname.EQ.'IKP')pnom(56)=pvalue
		if(pname.EQ.'TF')pnom(57)=pvalue
		if(pname.EQ.'QTF')pnom(58)=pvalue
		if(pname.EQ.'XTF')pnom(59)=pvalue
		if(pname.EQ.'VTF')pnom(60)=pvalue
		if(pname.EQ.'ITF')pnom(61)=pvalue
		if(pname.EQ.'TR')pnom(62)=pvalue
		if(pname.EQ.'TD')pnom(63)=pvalue
		if(pname.EQ.'KFN')pnom(64)=pvalue
		if(pname.EQ.'AFN')pnom(65)=pvalue
		if(pname.EQ.'BFN')pnom(66)=pvalue
		if(pname.EQ.'XRE')pnom(67)=pvalue
		if(pname.EQ.'XRBI')pnom(68)=pvalue
		if(pname.EQ.'XRCI')pnom(69)=pvalue
		if(pname.EQ.'XRS')pnom(70)=pvalue
		if(pname.EQ.'XVO')pnom(71)=pvalue
		if(pname.EQ.'EA')pnom(72)=pvalue
		if(pname.EQ.'EAIE')pnom(73)=pvalue
		if(pname.EQ.'EAIC')pnom(74)=pvalue
		if(pname.EQ.'EAIS')pnom(75)=pvalue
		if(pname.EQ.'EANE')pnom(76)=pvalue
		if(pname.EQ.'EANC')pnom(77)=pvalue
		if(pname.EQ.'EANS')pnom(78)=pvalue
		if(pname.EQ.'XIS')pnom(79)=pvalue
		if(pname.EQ.'XII')pnom(80)=pvalue
		if(pname.EQ.'XIN')pnom(81)=pvalue
		if(pname.EQ.'TNF')pnom(82)=pvalue
		if(pname.EQ.'TAVC')pnom(83)=pvalue
		if(pname.EQ.'RTH')pnom(84)=pvalue
		if(pname.EQ.'CTH')pnom(85)=pvalue
		if(pname.EQ.'VRT')pnom(86)=pvalue
		if(pname.EQ.'ART')pnom(87)=pvalue
		if(pname.EQ.'CCSO')pnom(88)=pvalue
		if(pname.EQ.'QBM')pnom(89)=pvalue
		if(pname.EQ.'NKF')pnom(90)=pvalue
		if(pname.EQ.'XIKF')pnom(91)=pvalue
		if(pname.EQ.'XRCX')pnom(92)=pvalue
		if(pname.EQ.'XRBX')pnom(93)=pvalue
		if(pname.EQ.'XRBP')pnom(94)=pvalue
		if(pname.EQ.'ISRR')pnom(95)=pvalue
		if(pname.EQ.'XISR')pnom(96)=pvalue
		if(pname.EQ.'DEAR')pnom(97)=pvalue
		if(pname.EQ.'EAP')pnom(98)=pvalue
		if(pname.EQ.'VBBE')pnom(99)=pvalue
		if(pname.EQ.'NBBE')pnom(100)=pvalue
		if(pname.EQ.'IBBE')pnom(101)=pvalue
		if(pname.EQ.'TVBBE1')pnom(102)=pvalue
		if(pname.EQ.'TVBBE2')pnom(103)=pvalue
		if(pname.EQ.'TNBBE')pnom(104)=pvalue
		if(pname.EQ.'EBBE')pnom(105)=pvalue
		if(pname.EQ.'DTEMP')pnom(106)=pvalue
		if(pname.EQ.'VERS')pnom(107)=pvalue
		if(pname.EQ.'VREV')pnom(108)=pvalue
		goto51
57		write(6,*) 'ERROR: cannot open file:',fname
		call exit(1)
58		continue
		close(31)
	endif
c
c	Check that resistances are not too small
c
	if(pnom(2).LT.1.0d-6)then
		write(6,*) 'ERROR: RCX too small'
		call exit(1)
	endif
	if(pnom(3).LT.1.0d-6)then
		write(6,*) 'ERROR: RCI too small'
		call exit(1)
	endif
	if(pnom(7).LT.1.0d-6)then
		write(6,*) 'ERROR: RBX too small'
		call exit(1)
	endif
	if(pnom(8).LT.1.0d-6)then
		write(6,*) 'ERROR: RBI too small'
		call exit(1)
	endif
	if(pnom(9).LT.1.0d-6)then
		write(6,*) 'ERROR: RE too small'
		call exit(1)
	endif
	if(pnom(10).LT.1.0d-6)then
		write(6,*) 'ERROR: RS too small'
		call exit(1)
	endif
	if(pnom(11).LT.1.0d-6)then
		write(6,*) 'ERROR: RBP too small'
		call exit(1)
	endif
	if(pnom(84).LT.1.0d-6)then
		write(6,*) 'ERROR: RTH too small'
		call exit(1)
	endif
c
c	Map parameters to desired temperature
c
	call vbic_3T_et_xf_t(p,pnom,TAMB)
c
c	Set up node initialization
c
	vtv=(273.15d0+p(1))*1.380662d-23/1.602189d-19
	vbinit=vtv*log(vtv/(1.414d0*p(12)))
	iinit=0
c
c	Now call the solver for each specified set of biases
c
100	continue
	if(doac.EQ.0)then
		read(5,*,end=200,err=199) (xx(i),i=1,3)
	else
		read(5,*,end=200,err=199) (xx(i),i=1,3),freq
	endif
101	continue
c
c	If needed, initialize node voltages
c
	xinit(1)=xx(1)
	xinit(2)=xx(1)
	xinit(3)=vbinit
	xinit(4)=vbinit
	xinit(5)=xx(3)
	xinit(6)=xx(1)
	xinit(7)=0.0d0
	xinit(8)=0.0d0
	xinit(9)=0.0d0
	xinit(10)=vbinit
	if(iinit.EQ.0)then
		do 102 i=1,10
			x(i)=xinit(i)
102		continue
		do 103 i=1,3
			xxlast(i)=xx(i)+1.0d0
103		continue
	endif
	dodc=0
	do 104 i=1,3
		if(xx(i).NE.xxlast(i))dodc=1
		xxlast(i)=xx(i)
104	continue
	if(dodc.NE.0)then
		call solver(p,xx,vmx,imx,x,nit,ierr,G,C,Frb,Fib,Frc,Fic)
		if(ierr.ne.0)then
c
c			If the solution fails, try continuation
c
			stm=1.0d0/8.0d0
			stp=stm
			val=0.0d0
			old=val
			do 118 i=1,10
				x(i)=xinit(i)
				xlast(i)=0.0d0
118			continue
			do 119 i=1,3
				cx(i)=val*xx(i)
119			continue
120			continue
			call solver(p,cx,vmx,imx,x,nit,ierr,G,C,Frb,Fib,Frc,Fic)
			if(ierr.EQ.0.AND.val.GE.1.0d0)goto140
			if(ierr.EQ.0)then
				old=val
				stp=stp*2.0d0
				if(stp.GT.stm)stp=stm
				val=val+stp
				if(val.GT.1.0d0)val=1.0d0
				do 131 i=1,10
					xlast(i)=x(i)
131				continue
				do 132 i=1,3
					cx(i)=xx(i)*val
132				continue
				goto120
			endif
			do 133 i=1,10
				x(i)=xlast(i)
133			continue
			stp=0.5*stp
			val=old+stp
			if(val.GT.1.0d0)val=1.0d0
			do 134 i=1,3
				cx(i)=xx(i)*val
134			continue
			if(stp.GE.2.0d-3)goto120
			write(6,*) 'ERROR: failed to converge'
			call exit(1)
140			continue
		endif
	else
		nit=0
		ierr=0
	endif
	iinit=1
	if(doac.EQ.0)then
		if(acc.EQ.0)then
			write(6,201) (vmx(i),i=1,4),nit
			write(6,201) (imx(i),i=1,3)
		else
			write(6,202) (vmx(i),i=1,4),nit
			write(6,202) (imx(i),i=1,3)
		endif
		goto100
	endif
c
c	AC solve is done based on MNA calculation, with
c	sequentially Vb=1+j0/Vc=0 and Vb=0/Vc=1+j0. For simplicity
c	real matrix solves only are done, not complex.
c	The Vb and Vc currents are introduced as unknowns.
c
	omega=2.0d0*3.141592654d0*freq
	do 150 i=1,11
		do 149 j=1,11
			Glu(i,j)=G(i,j)
			B(i,j)=omega*C(i,j)
149		continue
150	continue
	call decomp(Glu,nac,ipvtGlu,ierr)
	if(ierr.NE.0)then
		write(6,*) 'ERROR: failed AC solve'
		call exit(1)
	endif
	do 155 j=1,11
		do 153 i=1,11
			fr(i)=0.0d0
153		continue
		fr(j)=1.0d0
		call solve(Glu,nac,ipvtGlu,fr)
		do 154 i=1,11
			Ginv(i,j)=fr(i)
154		continue
155	continue
	do 160 i=1,11
		do 159 j=1,11
			BGinv(i,j)=0.0d0
			do 158 k=1,11
				BGinv(i,j)=BGinv(i,j)+B(i,k)*Ginv(k,j)
158			continue
159		continue
160	continue
	do 170 i=1,11
		do 169 j=1,11
			Geff(i,j)=G(i,j)
			do 168 k=1,11
				Geff(i,j)=Geff(i,j)+BGinv(i,k)*B(k,j)
168			continue
169		continue
170	continue
	call decomp(Geff,nac,ipvt,ierr)
	if(ierr.NE.0)then
		write(6,*) 'ERROR: failed AC solve'
		call exit(1)
	endif
	do 172 i=1,11
		fr(i)=Frb(i)
		do 171 j=1,11
			fr(i)=fr(i)+BGinv(i,j)*Fib(j)
171	continue
172	continue
	call solve(Geff,nac,ipvt,fr)
	y11r=fr(10)
	y21r=fr(11)
	do 174 i=1,11
		fi(i)=Fib(i)
		do 173 j=1,11
			fi(i)=fi(i)-B(i,j)*fr(j)
173	continue
174	continue
	call solve(Glu,nac,ipvtGlu,fi)
	y11i=fi(10)
	y21i=fi(11)
	do 182 i=1,11
		fr(i)=Frc(i)
		do 181 j=1,11
			fr(i)=fr(i)+BGinv(i,j)*Fic(j)
181	continue
182	continue
	call solve(Geff,nac,ipvt,fr)
	y12r=fr(10)
	y22r=fr(11)
	do 184 i=1,11
		fi(i)=Fic(i)
		do 183 j=1,11
			fi(i)=fi(i)-B(i,j)*fr(j)
183	continue
184	continue
	call solve(Glu,nac,ipvtGlu,fi)
	y12i=fi(10)
	y22i=fi(11)
	beta=sqrt((y21r*y21r+y21i*y21i)/(y11r*y11r+y11i*y11i))
	if(acc.EQ.0)then
		write(6,201) (vmx(i),i=1,4),nit
		write(6,201) (imx(i),i=1,3)
		write(6,201) freq
		write(6,201) beta
		write(6,201) y11r,y11i
		write(6,201) y12r,y12i
		write(6,201) y21r,y21i
		write(6,201) y22r,y22i
	else
		write(6,202) (vmx(i),i=1,4),nit
		write(6,202) (imx(i),i=1,3)
		write(6,202) freq
		write(6,202) beta
		write(6,202) y11r,y11i
		write(6,202) y12r,y12i
		write(6,202) y21r,y21i
		write(6,202) y22r,y22i
	endif
	goto100
199	continue
	write(6,*) 'ERROR: invalid bias specification'
	call exit(1)
200	continue
201	format(4(1x,1pe11.3),i4)
202	format(4(1x,1pe30.20),i4)
	call exit(0)
	end
c
c	Newton solver
c
	subroutine solver(p,xx,vmx,imx,x,nit,errflg,G,C,Gb,Cb,Gc,Cc)
	implicit double precision (a-z)
	integer n,i,j,k,ii,jj,itmax,nit,ipvt(10)
	integer errflg,nfail
	double precision p(108),xx(3),vmx(4),imx(3)
	double precision x(10),f(10),gdc(10,10)
	double precision G(11,11),C(11,11),xxlast(3)
	double precision Gb(11),Cb(11),Gc(11),Cc(11)
	itmax=200
	dvmax=0.2d0
	dtmax=10.0d0
	vatol=1.0d-8
	iatol=1.0d-12
	irtol=1.0e-3
	nit=0
	SCALE=1.0d0
	nfail=0
	n=10
c
c	No initialization, the system variables are stored
c	in the array x, and the previous solution is used
c	as the initial guess. This should be initialized
c	(to zero) on the first call. The order of unknowns is:
c	internal node     voltages
c	The formulation is nodal. Node order is:
c		cx ci bx bi ei bp xf1 xf2 [dt|b|dt b]
c
5	continue
10	continue

	Vbe=x(10)-xx(3)
	Vbc=x(10)-xx(1)
	Vbei=x(4)-x(5)
	Vbex=x(3)-x(5)
	Vbci=x(4)-x(2)
	Vbcx=x(4)-x(1)
	Vcei=x(2)-x(5)
	Veci=x(5)-x(2)
	Vbep=x(3)-x(6)
	Vre=xx(3)-x(5)
	Vrcx=xx(1)-x(1)
	Vrci=x(1)-x(2)
	Vrbx=x(10)-x(3)
	Vrbi=x(3)-x(4)
	Vrbp=x(6)-x(1)
	Vrth=x(9)
	Vith=x(9)
	Vizf=x(7)
	Vcxf=x(7)
	Vlxf=x(7)-x(8)
	Vrxf=x(8)
      call vbic_3T_et_xf_fj(p,Vrth,Vbei,Vbex,Vbci,Vbep,Vrxf,Vrcx
     +,Vbcx,Vrci,Vrbx,Vrbi,Vre,Vrbp,Vbe,Vbc,Vcei,Vcxf,Ibe,Ibe_Vrth
     +,Ibe_Vbei,Ibex,Ibex_Vrth,Ibex_Vbex,Itxf,Itxf_Vrxf,Itzr,Itzr_Vrth
     +,Itzr_Vbci,Itzr_Vbei,Ibc,Ibc_Vrth,Ibc_Vbci,Ibc_Vrxf,Ibc_Vbei,Ibep
     +,Ibep_Vrth,Ibep_Vbep,Ircx,Ircx_Vrcx,Ircx_Vrth,Irci,Irci_Vrci
     +,Irci_Vrth,Irci_Vbci,Irci_Vbcx,Irbx,Irbx_Vrbx,Irbx_Vrth,Irbi
     +,Irbi_Vrbi,Irbi_Vrth,Irbi_Vbei,Irbi_Vbci,Ire,Ire_Vre,Ire_Vrth,Irbp
     +,Irbp_Vrbp,Irbp_Vrth,Irbp_Vbep,Irbp_Vbci,Qbe,Qbe_Vrth,Qbe_Vbei
     +,Qbe_Vbci,Qbex,Qbex_Vrth,Qbex_Vbex,Qbc,Qbc_Vrth,Qbc_Vbci,Qbcx
     +,Qbcx_Vrth,Qbcx_Vbcx,Qbep,Qbep_Vrth,Qbep_Vbep,Qbep_Vbci,Qbeo
     +,Qbeo_Vbe,Qbco,Qbco_Vbc,Irth,Irth_Vrth,Ith,Ith_Vrth,Ith_Vbei
     +,Ith_Vbci,Ith_Vrxf,Ith_Vcei,Ith_Vbex,Ith_Vbep,Ith_Vrcx,Ith_Vrci
     +,Ith_Vbcx,Ith_Vrbx,Ith_Vrbi,Ith_Vre,Ith_Vrbp,Qcth,Qcth_Vrth,Ixzf
     +,Ixzf_Vrth,Ixzf_Vbei,Ixzf_Vbci,Ixxf,Ixxf_Vrxf,Qcxf,Qcxf_Vcxf,Flxf
     +,Flxf_Vrxf,SCALE)
	do 30 ii=1,10
		f(ii)=0.0d0
		do 25 jj=1,10
			gdc(ii,jj)=0.0d0
25		continue
30	continue
	iref=iatol/irtol
c
c	KCL at internal nodes
c
c
c	Stamp element: Ibe
c
	f(4)=f(4)-Ibe
	if(abs(Ibe).GT.iref)then
		iref=abs(Ibe)
	endif
	gdc(4,9)=gdc(4,9)-Ibe_Vrth
	gdc(4,4)=gdc(4,4)-Ibe_Vbei
	gdc(4,5)=gdc(4,5)+Ibe_Vbei
	f(5)=f(5)+Ibe
	gdc(5,9)=gdc(5,9)+Ibe_Vrth
	gdc(5,4)=gdc(5,4)+Ibe_Vbei
	gdc(5,5)=gdc(5,5)-Ibe_Vbei
c
c	Stamp element: Ibex
c
	f(3)=f(3)-Ibex
	if(abs(Ibex).GT.iref)then
		iref=abs(Ibex)
	endif
	gdc(3,9)=gdc(3,9)-Ibex_Vrth
	gdc(3,3)=gdc(3,3)-Ibex_Vbex
	gdc(3,5)=gdc(3,5)+Ibex_Vbex
	f(5)=f(5)+Ibex
	gdc(5,9)=gdc(5,9)+Ibex_Vrth
	gdc(5,3)=gdc(5,3)+Ibex_Vbex
	gdc(5,5)=gdc(5,5)-Ibex_Vbex
c
c	Stamp element: Itxf
c
	f(2)=f(2)-Itxf
	if(abs(Itxf).GT.iref)then
		iref=abs(Itxf)
	endif
	gdc(2,8)=gdc(2,8)-Itxf_Vrxf
	f(5)=f(5)+Itxf
	gdc(5,8)=gdc(5,8)+Itxf_Vrxf
c
c	Stamp element: Itzr
c
	f(5)=f(5)-Itzr
	if(abs(Itzr).GT.iref)then
		iref=abs(Itzr)
	endif
	gdc(5,9)=gdc(5,9)-Itzr_Vrth
	gdc(5,4)=gdc(5,4)-Itzr_Vbci
	gdc(5,2)=gdc(5,2)+Itzr_Vbci
	gdc(5,4)=gdc(5,4)-Itzr_Vbei
	gdc(5,5)=gdc(5,5)+Itzr_Vbei
	f(2)=f(2)+Itzr
	gdc(2,9)=gdc(2,9)+Itzr_Vrth
	gdc(2,4)=gdc(2,4)+Itzr_Vbci
	gdc(2,2)=gdc(2,2)-Itzr_Vbci
	gdc(2,4)=gdc(2,4)+Itzr_Vbei
	gdc(2,5)=gdc(2,5)-Itzr_Vbei
c
c	Stamp element: Ibc
c
	f(4)=f(4)-Ibc
	if(abs(Ibc).GT.iref)then
		iref=abs(Ibc)
	endif
	gdc(4,9)=gdc(4,9)-Ibc_Vrth
	gdc(4,4)=gdc(4,4)-Ibc_Vbci
	gdc(4,2)=gdc(4,2)+Ibc_Vbci
	gdc(4,8)=gdc(4,8)-Ibc_Vrxf
	gdc(4,4)=gdc(4,4)-Ibc_Vbei
	gdc(4,5)=gdc(4,5)+Ibc_Vbei
	f(2)=f(2)+Ibc
	gdc(2,9)=gdc(2,9)+Ibc_Vrth
	gdc(2,4)=gdc(2,4)+Ibc_Vbci
	gdc(2,2)=gdc(2,2)-Ibc_Vbci
	gdc(2,8)=gdc(2,8)+Ibc_Vrxf
	gdc(2,4)=gdc(2,4)+Ibc_Vbei
	gdc(2,5)=gdc(2,5)-Ibc_Vbei
c
c	Stamp element: Ibep
c
	f(3)=f(3)-Ibep
	if(abs(Ibep).GT.iref)then
		iref=abs(Ibep)
	endif
	gdc(3,9)=gdc(3,9)-Ibep_Vrth
	gdc(3,3)=gdc(3,3)-Ibep_Vbep
	gdc(3,6)=gdc(3,6)+Ibep_Vbep
	f(6)=f(6)+Ibep
	gdc(6,9)=gdc(6,9)+Ibep_Vrth
	gdc(6,3)=gdc(6,3)+Ibep_Vbep
	gdc(6,6)=gdc(6,6)-Ibep_Vbep
c
c	Stamp element: Ircx
c
	f(1)=f(1)+Ircx
	gdc(1,1)=gdc(1,1)-Ircx_Vrcx
	gdc(1,9)=gdc(1,9)+Ircx_Vrth
c
c	Stamp element: Irci
c
	f(1)=f(1)-Irci
	if(abs(Irci).GT.iref)then
		iref=abs(Irci)
	endif
	gdc(1,1)=gdc(1,1)-Irci_Vrci
	gdc(1,2)=gdc(1,2)+Irci_Vrci
	gdc(1,9)=gdc(1,9)-Irci_Vrth
	gdc(1,4)=gdc(1,4)-Irci_Vbci
	gdc(1,2)=gdc(1,2)+Irci_Vbci
	gdc(1,4)=gdc(1,4)-Irci_Vbcx
	gdc(1,1)=gdc(1,1)+Irci_Vbcx
	f(2)=f(2)+Irci
	gdc(2,1)=gdc(2,1)+Irci_Vrci
	gdc(2,2)=gdc(2,2)-Irci_Vrci
	gdc(2,9)=gdc(2,9)+Irci_Vrth
	gdc(2,4)=gdc(2,4)+Irci_Vbci
	gdc(2,2)=gdc(2,2)-Irci_Vbci
	gdc(2,4)=gdc(2,4)+Irci_Vbcx
	gdc(2,1)=gdc(2,1)-Irci_Vbcx
c
c	Stamp element: Irbx
c
	f(10)=f(10)-Irbx
	if(abs(Irbx).GT.iref)then
		iref=abs(Irbx)
	endif
	gdc(10,10)=gdc(10,10)-Irbx_Vrbx
	gdc(10,3)=gdc(10,3)+Irbx_Vrbx
	gdc(10,9)=gdc(10,9)-Irbx_Vrth
	f(3)=f(3)+Irbx
	gdc(3,10)=gdc(3,10)+Irbx_Vrbx
	gdc(3,3)=gdc(3,3)-Irbx_Vrbx
	gdc(3,9)=gdc(3,9)+Irbx_Vrth
c
c	Stamp element: Irbi
c
	f(3)=f(3)-Irbi
	if(abs(Irbi).GT.iref)then
		iref=abs(Irbi)
	endif
	gdc(3,3)=gdc(3,3)-Irbi_Vrbi
	gdc(3,4)=gdc(3,4)+Irbi_Vrbi
	gdc(3,9)=gdc(3,9)-Irbi_Vrth
	gdc(3,4)=gdc(3,4)-Irbi_Vbei
	gdc(3,5)=gdc(3,5)+Irbi_Vbei
	gdc(3,4)=gdc(3,4)-Irbi_Vbci
	gdc(3,2)=gdc(3,2)+Irbi_Vbci
	f(4)=f(4)+Irbi
	gdc(4,3)=gdc(4,3)+Irbi_Vrbi
	gdc(4,4)=gdc(4,4)-Irbi_Vrbi
	gdc(4,9)=gdc(4,9)+Irbi_Vrth
	gdc(4,4)=gdc(4,4)+Irbi_Vbei
	gdc(4,5)=gdc(4,5)-Irbi_Vbei
	gdc(4,4)=gdc(4,4)+Irbi_Vbci
	gdc(4,2)=gdc(4,2)-Irbi_Vbci
c
c	Stamp element: Ire
c
	f(5)=f(5)+Ire
	gdc(5,5)=gdc(5,5)-Ire_Vre
	gdc(5,9)=gdc(5,9)+Ire_Vrth
c
c	Stamp element: Irbp
c
	f(6)=f(6)-Irbp
	if(abs(Irbp).GT.iref)then
		iref=abs(Irbp)
	endif
	gdc(6,6)=gdc(6,6)-Irbp_Vrbp
	gdc(6,1)=gdc(6,1)+Irbp_Vrbp
	gdc(6,9)=gdc(6,9)-Irbp_Vrth
	gdc(6,3)=gdc(6,3)-Irbp_Vbep
	gdc(6,6)=gdc(6,6)+Irbp_Vbep
	gdc(6,4)=gdc(6,4)-Irbp_Vbci
	gdc(6,2)=gdc(6,2)+Irbp_Vbci
	f(1)=f(1)+Irbp
	gdc(1,6)=gdc(1,6)+Irbp_Vrbp
	gdc(1,1)=gdc(1,1)-Irbp_Vrbp
	gdc(1,9)=gdc(1,9)+Irbp_Vrth
	gdc(1,3)=gdc(1,3)+Irbp_Vbep
	gdc(1,6)=gdc(1,6)-Irbp_Vbep
	gdc(1,4)=gdc(1,4)+Irbp_Vbci
	gdc(1,2)=gdc(1,2)-Irbp_Vbci
c
c	Stamp element: Irth
c
	f(9)=f(9)-Irth
	gdc(9,9)=gdc(9,9)-Irth_Vrth
c
c	Stamp element: Ith
c
	f(9)=f(9)-Ith
	gdc(9,9)=gdc(9,9)-Ith_Vrth
	gdc(9,4)=gdc(9,4)-Ith_Vbei
	gdc(9,5)=gdc(9,5)+Ith_Vbei
	gdc(9,4)=gdc(9,4)-Ith_Vbci
	gdc(9,2)=gdc(9,2)+Ith_Vbci
	gdc(9,8)=gdc(9,8)-Ith_Vrxf
	gdc(9,2)=gdc(9,2)-Ith_Vcei
	gdc(9,5)=gdc(9,5)+Ith_Vcei
	gdc(9,3)=gdc(9,3)-Ith_Vbex
	gdc(9,5)=gdc(9,5)+Ith_Vbex
	gdc(9,3)=gdc(9,3)-Ith_Vbep
	gdc(9,6)=gdc(9,6)+Ith_Vbep
	gdc(9,1)=gdc(9,1)+Ith_Vrcx
	gdc(9,1)=gdc(9,1)-Ith_Vrci
	gdc(9,2)=gdc(9,2)+Ith_Vrci
	gdc(9,4)=gdc(9,4)-Ith_Vbcx
	gdc(9,1)=gdc(9,1)+Ith_Vbcx
	gdc(9,10)=gdc(9,10)-Ith_Vrbx
	gdc(9,3)=gdc(9,3)+Ith_Vrbx
	gdc(9,3)=gdc(9,3)-Ith_Vrbi
	gdc(9,4)=gdc(9,4)+Ith_Vrbi
	gdc(9,5)=gdc(9,5)+Ith_Vre
	gdc(9,6)=gdc(9,6)-Ith_Vrbp
	gdc(9,1)=gdc(9,1)+Ith_Vrbp
c
c	Stamp element: Ixzf
c
	f(7)=f(7)-Ixzf
	if(abs(Ixzf).GT.iref)then
		iref=abs(Ixzf)
	endif
	gdc(7,9)=gdc(7,9)-Ixzf_Vrth
	gdc(7,4)=gdc(7,4)-Ixzf_Vbei
	gdc(7,5)=gdc(7,5)+Ixzf_Vbei
	gdc(7,4)=gdc(7,4)-Ixzf_Vbci
	gdc(7,2)=gdc(7,2)+Ixzf_Vbci
c
c	Stamp element: Ixxf
c
	f(8)=f(8)-Ixxf
	if(abs(Ixxf).GT.iref)then
		iref=abs(Ixxf)
	endif
	gdc(8,8)=gdc(8,8)-Ixxf_Vrxf
c
c	Stamp element: ib
c
	f(10)=f(10)+xx(2)
c
c	Stamp xf elements (NOT MNA, reduced form)
c
	f(7)=f(7)-x(8)
	gdc(7,8)=gdc(7,8)-1.0d0
	f(8)=f(8)+x(7)
	gdc(8,7)=gdc(8,7)+1.0d0
c
c	Solve for update vector
c
	ianorm=0.0d0
	do 80 i=1,10
		if(abs(f(i)).GT.ianorm)ianorm=abs(f(i))
80	continue
	irnorm=0.0d0
	do 81 i=1,10
		relerr=abs(f(i))/iref
		if(relerr.GT.irnorm)irnorm=relerr
81	continue
	call decomp(gdc,n,ipvt,errflg)
	if(errflg.NE.0)then
		if(nfail.EQ.0)then
c
c		On first failure, try starting from zero
c
			do 88 i=1,10
				x(i)=0.0d0
88			continue
			nit=0
			nfail=nfail+1
			goto5
		endif
		return
	endif
	call solve(gdc,n,ipvt,f)
	vnorm=0.0d0
	vscale=1.0d0
	do 90 i=1,10
		if(i.EQ.7)goto90
		if(i.EQ.8)goto90
		if(abs(f(i)).GT.vnorm)vnorm=abs(f(i))
		if(i.EQ.9)then
			if(abs(f(i)/dtmax).GT.1.0d0/vscale)vscale=abs(dtmax/f(i))
		else
			if(abs(f(i)/dvmax).GT.1.0d0/vscale)vscale=abs(dvmax/f(i))
		endif
90	continue
	do 100 i=1,10
		x(i)=x(i)-(f(i)*vscale)
100	continue
	nit=nit+1
	if(nit.LT.itmax.AND.(ianorm.GT.iatol.AND.irnorm.GT.irtol))goto10
	if(nit.LT.itmax.AND.(vnorm.GT.vatol))goto10
	if(nit.LT.3)goto10
	if(nit.EQ.itmax)errflg=1
c
c	Copy terminal solution into arrays to return
c
	do 200 i=1,3
		vmx(i)=xx(i)
		imx(i)=0.0d0
200	continue
	vmx(4)=x(9)
	vmx(2)=x(10)
	imx(1)=imx(1)+Ircx
	imx(2)=imx(2)+Irbx
	imx(3)=imx(3)+Ire
c
c	Finally generate complete AC matrices, set up for y-parameter calculation
c
	do 300 ii=1,11
		Gb(ii)=0.0d0
		Cb(ii)=0.0d0
		Gc(ii)=0.0d0
		Cc(ii)=0.0d0
		do 295 jj=1,11
			G(ii,jj)=0.0d0
			C(ii,jj)=0.0d0
295		continue
300	continue
c
c	Stamp element: Ibe
c
	G(4,9)=G(4,9)-Ibe_Vrth
	G(4,4)=G(4,4)-Ibe_Vbei
	G(4,5)=G(4,5)+Ibe_Vbei
	G(5,9)=G(5,9)+Ibe_Vrth
	G(5,4)=G(5,4)+Ibe_Vbei
	G(5,5)=G(5,5)-Ibe_Vbei
c
c	Stamp element: Ibex
c
	G(3,9)=G(3,9)-Ibex_Vrth
	G(3,3)=G(3,3)-Ibex_Vbex
	G(3,5)=G(3,5)+Ibex_Vbex
	G(5,9)=G(5,9)+Ibex_Vrth
	G(5,3)=G(5,3)+Ibex_Vbex
	G(5,5)=G(5,5)-Ibex_Vbex
c
c	Stamp element: Itxf
c
	G(2,8)=G(2,8)-Itxf_Vrxf
	G(5,8)=G(5,8)+Itxf_Vrxf
c
c	Stamp element: Itzr
c
	G(5,9)=G(5,9)-Itzr_Vrth
	G(5,4)=G(5,4)-Itzr_Vbci
	G(5,2)=G(5,2)+Itzr_Vbci
	G(5,4)=G(5,4)-Itzr_Vbei
	G(5,5)=G(5,5)+Itzr_Vbei
	G(2,9)=G(2,9)+Itzr_Vrth
	G(2,4)=G(2,4)+Itzr_Vbci
	G(2,2)=G(2,2)-Itzr_Vbci
	G(2,4)=G(2,4)+Itzr_Vbei
	G(2,5)=G(2,5)-Itzr_Vbei
c
c	Stamp element: Ibc
c
	G(4,9)=G(4,9)-Ibc_Vrth
	G(4,4)=G(4,4)-Ibc_Vbci
	G(4,2)=G(4,2)+Ibc_Vbci
	G(4,8)=G(4,8)-Ibc_Vrxf
	G(4,4)=G(4,4)-Ibc_Vbei
	G(4,5)=G(4,5)+Ibc_Vbei
	G(2,9)=G(2,9)+Ibc_Vrth
	G(2,4)=G(2,4)+Ibc_Vbci
	G(2,2)=G(2,2)-Ibc_Vbci
	G(2,8)=G(2,8)+Ibc_Vrxf
	G(2,4)=G(2,4)+Ibc_Vbei
	G(2,5)=G(2,5)-Ibc_Vbei
c
c	Stamp element: Ibep
c
	G(3,9)=G(3,9)-Ibep_Vrth
	G(3,3)=G(3,3)-Ibep_Vbep
	G(3,6)=G(3,6)+Ibep_Vbep
	G(6,9)=G(6,9)+Ibep_Vrth
	G(6,3)=G(6,3)+Ibep_Vbep
	G(6,6)=G(6,6)-Ibep_Vbep
c
c	Stamp element: Ircx
c
	Gc(11)=Gc(11)+Ircx_Vrcx
	G(11,1)=G(11,1)+Ircx_Vrcx
	G(11,9)=G(11,9)-Ircx_Vrth
	Gc(1)=Gc(1)-Ircx_Vrcx
	G(1,1)=G(1,1)-Ircx_Vrcx
	G(1,9)=G(1,9)+Ircx_Vrth
c
c	Stamp element: Irci
c
	G(1,1)=G(1,1)-Irci_Vrci
	G(1,2)=G(1,2)+Irci_Vrci
	G(1,9)=G(1,9)-Irci_Vrth
	G(1,4)=G(1,4)-Irci_Vbci
	G(1,2)=G(1,2)+Irci_Vbci
	G(1,4)=G(1,4)-Irci_Vbcx
	G(1,1)=G(1,1)+Irci_Vbcx
	G(2,1)=G(2,1)+Irci_Vrci
	G(2,2)=G(2,2)-Irci_Vrci
	G(2,9)=G(2,9)+Irci_Vrth
	G(2,4)=G(2,4)+Irci_Vbci
	G(2,2)=G(2,2)-Irci_Vbci
	G(2,4)=G(2,4)+Irci_Vbcx
	G(2,1)=G(2,1)-Irci_Vbcx
c
c	Stamp element: Irbx
c
	Gb(10)=Gb(10)+Irbx_Vrbx
	G(10,3)=G(10,3)+Irbx_Vrbx
	G(10,9)=G(10,9)-Irbx_Vrth
	Gb(3)=Gb(3)-Irbx_Vrbx
	G(3,3)=G(3,3)-Irbx_Vrbx
	G(3,9)=G(3,9)+Irbx_Vrth
c
c	Stamp element: Irbi
c
	G(3,3)=G(3,3)-Irbi_Vrbi
	G(3,4)=G(3,4)+Irbi_Vrbi
	G(3,9)=G(3,9)-Irbi_Vrth
	G(3,4)=G(3,4)-Irbi_Vbei
	G(3,5)=G(3,5)+Irbi_Vbei
	G(3,4)=G(3,4)-Irbi_Vbci
	G(3,2)=G(3,2)+Irbi_Vbci
	G(4,3)=G(4,3)+Irbi_Vrbi
	G(4,4)=G(4,4)-Irbi_Vrbi
	G(4,9)=G(4,9)+Irbi_Vrth
	G(4,4)=G(4,4)+Irbi_Vbei
	G(4,5)=G(4,5)-Irbi_Vbei
	G(4,4)=G(4,4)+Irbi_Vbci
	G(4,2)=G(4,2)-Irbi_Vbci
c
c	Stamp element: Ire
c
	G(5,5)=G(5,5)-Ire_Vre
	G(5,9)=G(5,9)+Ire_Vrth
c
c	Stamp element: Irbp
c
	G(6,6)=G(6,6)-Irbp_Vrbp
	G(6,1)=G(6,1)+Irbp_Vrbp
	G(6,9)=G(6,9)-Irbp_Vrth
	G(6,3)=G(6,3)-Irbp_Vbep
	G(6,6)=G(6,6)+Irbp_Vbep
	G(6,4)=G(6,4)-Irbp_Vbci
	G(6,2)=G(6,2)+Irbp_Vbci
	G(1,6)=G(1,6)+Irbp_Vrbp
	G(1,1)=G(1,1)-Irbp_Vrbp
	G(1,9)=G(1,9)+Irbp_Vrth
	G(1,3)=G(1,3)+Irbp_Vbep
	G(1,6)=G(1,6)-Irbp_Vbep
	G(1,4)=G(1,4)+Irbp_Vbci
	G(1,2)=G(1,2)-Irbp_Vbci
c
c	Stamp element: Irth
c
	G(9,9)=G(9,9)-Irth_Vrth
c
c	Stamp element: Ith
c
	G(9,9)=G(9,9)-Ith_Vrth
	G(9,4)=G(9,4)-Ith_Vbei
	G(9,5)=G(9,5)+Ith_Vbei
	G(9,4)=G(9,4)-Ith_Vbci
	G(9,2)=G(9,2)+Ith_Vbci
	G(9,8)=G(9,8)-Ith_Vrxf
	G(9,2)=G(9,2)-Ith_Vcei
	G(9,5)=G(9,5)+Ith_Vcei
	G(9,3)=G(9,3)-Ith_Vbex
	G(9,5)=G(9,5)+Ith_Vbex
	G(9,3)=G(9,3)-Ith_Vbep
	G(9,6)=G(9,6)+Ith_Vbep
	G(9,1)=G(9,1)+Ith_Vrcx
	G(9,1)=G(9,1)-Ith_Vrci
	G(9,2)=G(9,2)+Ith_Vrci
	G(9,4)=G(9,4)-Ith_Vbcx
	G(9,1)=G(9,1)+Ith_Vbcx
	G(9,3)=G(9,3)+Ith_Vrbx
	G(9,3)=G(9,3)-Ith_Vrbi
	G(9,4)=G(9,4)+Ith_Vrbi
	G(9,5)=G(9,5)+Ith_Vre
	G(9,6)=G(9,6)-Ith_Vrbp
	G(9,1)=G(9,1)+Ith_Vrbp
c
c	Stamp element: Ixzf
c
	G(7,9)=G(7,9)-Ixzf_Vrth
	G(7,4)=G(7,4)-Ixzf_Vbei
	G(7,5)=G(7,5)+Ixzf_Vbei
	G(7,4)=G(7,4)-Ixzf_Vbci
	G(7,2)=G(7,2)+Ixzf_Vbci
c
c	Stamp element: Ixxf
c
	G(8,8)=G(8,8)-Ixxf_Vrxf
c
c	Stamp element: Qbe
c
	C(4,9)=C(4,9)-Qbe_Vrth
	C(4,4)=C(4,4)-Qbe_Vbei
	C(4,5)=C(4,5)+Qbe_Vbei
	C(4,4)=C(4,4)-Qbe_Vbci
	C(4,2)=C(4,2)+Qbe_Vbci
	C(5,9)=C(5,9)+Qbe_Vrth
	C(5,4)=C(5,4)+Qbe_Vbei
	C(5,5)=C(5,5)-Qbe_Vbei
	C(5,4)=C(5,4)+Qbe_Vbci
	C(5,2)=C(5,2)-Qbe_Vbci
c
c	Stamp element: Qbex
c
	C(3,9)=C(3,9)-Qbex_Vrth
	C(3,3)=C(3,3)-Qbex_Vbex
	C(3,5)=C(3,5)+Qbex_Vbex
	C(5,9)=C(5,9)+Qbex_Vrth
	C(5,3)=C(5,3)+Qbex_Vbex
	C(5,5)=C(5,5)-Qbex_Vbex
c
c	Stamp element: Qbc
c
	C(4,9)=C(4,9)-Qbc_Vrth
	C(4,4)=C(4,4)-Qbc_Vbci
	C(4,2)=C(4,2)+Qbc_Vbci
	C(2,9)=C(2,9)+Qbc_Vrth
	C(2,4)=C(2,4)+Qbc_Vbci
	C(2,2)=C(2,2)-Qbc_Vbci
c
c	Stamp element: Qbcx
c
	C(4,9)=C(4,9)-Qbcx_Vrth
	C(4,4)=C(4,4)-Qbcx_Vbcx
	C(4,1)=C(4,1)+Qbcx_Vbcx
	C(1,9)=C(1,9)+Qbcx_Vrth
	C(1,4)=C(1,4)+Qbcx_Vbcx
	C(1,1)=C(1,1)-Qbcx_Vbcx
c
c	Stamp element: Qbep
c
	C(3,9)=C(3,9)-Qbep_Vrth
	C(3,3)=C(3,3)-Qbep_Vbep
	C(3,6)=C(3,6)+Qbep_Vbep
	C(3,4)=C(3,4)-Qbep_Vbci
	C(3,2)=C(3,2)+Qbep_Vbci
	C(6,9)=C(6,9)+Qbep_Vrth
	C(6,3)=C(6,3)+Qbep_Vbep
	C(6,6)=C(6,6)-Qbep_Vbep
	C(6,4)=C(6,4)+Qbep_Vbci
	C(6,2)=C(6,2)-Qbep_Vbci
c
c	Stamp element: Qbeo
c
	Cb(10)=Cb(10)+Qbeo_Vbe
c
c	Stamp element: Qbco
c
	Cb(10)=Cb(10)+Qbco_Vbc
	Cc(10)=Cc(10)-Qbco_Vbc
	Cb(11)=Cb(11)-Qbco_Vbc
	Cb(11)=Cb(11)+Qbco_Vbc
c
c	Stamp element: Qcth
c
	C(9,9)=C(9,9)-Qcth_Vrth
c
c	Stamp element: Qcxf
c
	C(7,7)=C(7,7)-Qcxf_Vcxf
c
c	Stamp xf elements (NOT MNA, reduced form)
c
	G(7,8)=G(7,8)-1.0d0
	G(8,7)=G(8,7)+1.0d0
	C(8,8)=C(8,8)-Flxf_Vrxf
	G(10,10)=G(10,10)+1.0d0
	G(11,11)=G(11,11)+1.0d0
	return
	end
c
c	FMM forward decomposition and back substitution routines
c
	subroutine decomp(a,n,ipvt,errflg)
	integer n,ipvt(n),errflg
	integer nm1,k,kp1,m,i,j
	double precision a(n,n)
	double precision t
	errflg=0
	ipvt(n)=1
	nm1 = n-1
	do 35 k=1,nm1
		kp1=k+1
		m=k
		do 15 i=kp1,n
			if(abs(a(i,k)).gt.abs(a(m,k))) m=i
15		continue
		ipvt(k)=m
		if(m.ne.k) ipvt(n)=-ipvt(n)
		t=a(m,k)
		a(m,k)=a(k,k)
		a(k,k)=t
		if(t .eq. 0.0) then
			errflg=1
			return
		endif
		do 20 i=kp1,n
			a(i,k) = -a(i,k)/t
20		continue
		do 30 j=kp1,n
			t=a(m,j)
			a(m,j)=a(k,j)
			a(k,j)=t
			do 25 i=kp1,n
				a(i,j)=a(i,j)+a(i,k)*t
25			continue
30		continue
35	continue
	return
	end
	subroutine solve(a,n,ipvt,b)
	integer n,ipvt(n)
	double precision a(n,n),b(n)
	integer nm1,k,kp1,m,i,kb,km1
	double precision t
	nm1=n-1
	do 20 k=1,nm1
		kp1=k+1
		m=ipvt(k)
		t=b(m)
		b(m)=b(k)
		b(k)=t
		do 10 i=kp1,n
			b(i)=b(i)+a(i,k)*t
10		continue
20	continue
	do 40 kb=1,nm1
		km1=n-kb
		k=km1+1
		b(k)=b(k)/a(k,k)
		t=-b(k)
		do 30 i=1,km1
			b(i)=b(i)+a(i,k)*t
30		continue
40	continue
	b(1)=b(1)/a(1,1)
	return
	end
