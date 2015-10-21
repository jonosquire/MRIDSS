// Automatically generated temporary variables - class definition
dcmplxVecM T2TdzP2TiLap2Tkxb, T2TdzTiLap2Tky, T2TdzTOm, TdzTiLap2Tkxb, TdzTkxb, TdzTSPLUSTm2TdzTOm, TiLap2TkxbP2Tky, TiLap2Tky, Tm2TdzTiLap2TkxbP2Tky, Tm2TiLap2TkxbTkyP2, TmdzTiLap2Tkxb, TmdzTiLap2TkxbP3, TmiLap2TkxbP2Tky, TmkxbTkyTqPLUSTmkxbTkyTS;

dcmplx Tkxb, Tky, Tmkxb, Tmky;

dcmplxMat Ctmp_1_, Ctmp_2_, Ctmp_3_, Ctmp_4_, Ctmp_5_, Ctmp_6_;

dcmplxMat C11_, C12_, C21_, C22_;






// Automatically generated temporary variables - class constructor
T2TdzP2TiLap2Tkxb = dcmplxVecM(NZ_);
T2TdzTiLap2Tky = dcmplxVecM(NZ_);
T2TdzTOm = dcmplxVecM(NZ_);
TdzTiLap2Tkxb = dcmplxVecM(NZ_);
TdzTkxb = dcmplxVecM(NZ_);
TdzTSPLUSTm2TdzTOm = dcmplxVecM(NZ_);
TiLap2TkxbP2Tky = dcmplxVecM(NZ_);
TiLap2Tky = dcmplxVecM(NZ_);
Tm2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
Tm2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TmdzTiLap2Tkxb = dcmplxVecM(NZ_);
TmdzTiLap2TkxbP3 = dcmplxVecM(NZ_);
TmiLap2TkxbP2Tky = dcmplxVecM(NZ_);
TmkxbTkyTqPLUSTmkxbTkyTS = dcmplxVecM(NZ_);

Ctmp_1_ = dcmplxMat(NZ_,NZ_);
Ctmp_2_ = dcmplxMat(NZ_,NZ_);
Ctmp_3_ = dcmplxMat(NZ_,NZ_);
Ctmp_4_ = dcmplxMat(NZ_,NZ_);
Ctmp_5_ = dcmplxMat(NZ_,NZ_);
Ctmp_6_ = dcmplxMat(NZ_,NZ_);

C11_ = dcmplxMat(NZ_,NZ_);
C12_ = dcmplxMat(NZ_,NZ_);
C21_ = dcmplxMat(NZ_,NZ_);
C22_ = dcmplxMat(NZ_,NZ_);








//Assign automatically generated variables in equations 
// Submatrix variable definition (automatic)
C11_ = Ckl_in[i].block( 0, 0, NZ_, NZ_);
C12_ = Ckl_in[i].block( 0, NZ_, NZ_, NZ_);
C21_ = Ckl_in[i].block( NZ_, 0, NZ_, NZ_);
C22_ = Ckl_in[i].block( NZ_, NZ_, NZ_, NZ_);

// Scalar variable definition (automatic)
Tkxb = (kxctmp_)*fft2Dfac_;
Tky = (kyctmp_)*fft2Dfac_;
Tmkxb = ((-1.*kxctmp_))*fft2Dfac_;
Tmky = ((-1.*kyctmp_))*fft2Dfac_;

// Vector variable definition (automatic)
T2TdzP2TiLap2Tkxb = ((2.*kxctmp_)*(ilap2tmp_*kz2_).matrix())*fft2Dfac_;
T2TdzTiLap2Tky = ((2.*kyctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
T2TdzTOm = ((2.*Omega_)*kz_.matrix())*fft2Dfac_;
TdzTiLap2Tkxb = ((ilap2tmp_*kz_).matrix()*kxctmp_)*fft2Dfac_;
TdzTkxb = (kxctmp_*kz_.matrix())*fft2Dfac_;
TdzTSPLUSTm2TdzTOm = ((-2.*Omega_)*kz_.matrix() + kz_.matrix()*q_)*fft2Dfac_;
TiLap2TkxbP2Tky = (ilap2tmp_.matrix()*(kyctmp_*pow(kxctmp_,2)))*fft2Dfac_;
TiLap2Tky = (ilap2tmp_.matrix()*kyctmp_)*fft2Dfac_;
Tm2TdzTiLap2TkxbP2Tky = ((-2.*(kyctmp_*pow(kxctmp_,2)))*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
Tm2TiLap2TkxbTkyP2 = ((-2.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
TmdzTiLap2Tkxb = ((-1.*kxctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
TmdzTiLap2TkxbP3 = ((-1.*pow(kxctmp_,3))*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
TmiLap2TkxbP2Tky = ((-1.*(kyctmp_*pow(kxctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
TmkxbTkyTqPLUSTmkxbTkyTS = ((-2.*(kxctmp_*kyctmp_))*q_)*fft2Dfac_;







////////////////////////////////////////////////////////
///       AUTOMATICALLY GENERATED EQUATIONS         /////
///       see GenerateC++Equations.nb in MMA        /////

Ctmp_1_ = Tmkxb*C11_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Ux_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tmky*C11_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = Uy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C11_ + T2TdzTiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzUx_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = dzdzdzUx_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C11_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzUy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTkxb.asDiagonal()*C11_ + Tky*C21_ + TmdzTiLap2TkxbP3.asDiagonal()*C11_ + TmiLap2TkxbP2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
fft_.for_2DFull( Ctmp_3_ );
Ctmp_6_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_3_ + (T2TdzTOm/fft2Dfac_).asDiagonal()*C21_ + (TmkxbTkyTqPLUSTmkxbTkyTS/fft2Dfac_).asDiagonal()*C11_);
dealias(Ctmp_6_);
Ckl_out[i].block( 0, 0, NZ_, NZ_) = Ctmp_6_;



Ctmp_1_ = Tmkxb*C12_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Ux_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tmky*C12_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = Uy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C12_ + T2TdzTiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzUx_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = dzdzdzUx_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C12_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzUy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTkxb.asDiagonal()*C12_ + Tky*C22_ + TmdzTiLap2TkxbP3.asDiagonal()*C12_ + TmiLap2TkxbP2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
fft_.for_2DFull( Ctmp_3_ );
Ctmp_6_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_3_ + (T2TdzTOm/fft2Dfac_).asDiagonal()*C22_ + (TmkxbTkyTqPLUSTmkxbTkyTS/fft2Dfac_).asDiagonal()*C12_);
dealias(Ctmp_6_);
Ckl_out[i].block( 0, NZ_, NZ_, NZ_) = Ctmp_6_;



Ctmp_1_ = Tmkxb*C21_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Ux_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tmky*C21_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = Uy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzUy_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = TiLap2TkxbP2Tky.asDiagonal()*C11_ + TmdzTiLap2Tkxb.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = dzUx_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = Tkxb*C11_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzUy_.asDiagonal()*Ctmp_5_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_5_ = Ctmp_1_ + (TdzTSPLUSTm2TdzTOm/fft2Dfac_).asDiagonal()*C11_;
dealias(Ctmp_5_);
Ckl_out[i].block( NZ_, 0, NZ_, NZ_) = Ctmp_5_;



Ctmp_1_ = Tmkxb*C22_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Ux_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tmky*C22_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = Uy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzUy_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = TiLap2TkxbP2Tky.asDiagonal()*C12_ + TmdzTiLap2Tkxb.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = dzUx_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = Tkxb*C12_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzUy_.asDiagonal()*Ctmp_5_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_5_ = Ctmp_1_ + (TdzTSPLUSTm2TdzTOm/fft2Dfac_).asDiagonal()*C12_;
dealias(Ctmp_5_);
Ckl_out[i].block( NZ_, NZ_, NZ_, NZ_) = Ctmp_5_;
///     AUTOMATICALLY GENERATED EQUATIONS - END       ///
////////////////////////////////////////////////////////
