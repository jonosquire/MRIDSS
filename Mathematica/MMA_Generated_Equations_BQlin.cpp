// Automatically generated temporary variables - class definition
dcmplxVecM T2Tdz, T2TdzTiLap2TkxbP2Tky, T2TdzTiLap2Tky, T2TiLap2TkxbTkyP2, TB0zTdz, TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2, TdzTiLap2Tkxb, TdzTqPLUSTm2Tdz, TiLap2Tky, TmdzTiLap2Tkxb, TmdzTq, TmiLap2Tky;

dcmplx Tky, Tm2TkxbTkyTq, Tmkxb;

dcmplxMat Ctmp_1_, Ctmp_2_, Ctmp_3_;

dcmplxMat C11_, C12_, C13_, C14_, C21_, C22_, C23_, C24_, C31_, C32_, C33_, C34_, C41_, C42_, C43_, C44_;






// Automatically generated temporary variables - class constructor
T2Tdz = dcmplxVecM(NZ_);
T2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
T2TdzTiLap2Tky = dcmplxVecM(NZ_);
T2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TB0zTdz = dcmplxVecM(NZ_);
TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TdzTiLap2Tkxb = dcmplxVecM(NZ_);
TdzTqPLUSTm2Tdz = dcmplxVecM(NZ_);
TiLap2Tky = dcmplxVecM(NZ_);
TmdzTiLap2Tkxb = dcmplxVecM(NZ_);
TmdzTq = dcmplxVecM(NZ_);
TmiLap2Tky = dcmplxVecM(NZ_);

Ctmp_1_ = dcmplxMat(NZ_,NZ_);
Ctmp_2_ = dcmplxMat(NZ_,NZ_);
Ctmp_3_ = dcmplxMat(NZ_,NZ_);

C11_ = dcmplxMat(NZ_,NZ_);
C12_ = dcmplxMat(NZ_,NZ_);
C13_ = dcmplxMat(NZ_,NZ_);
C14_ = dcmplxMat(NZ_,NZ_);
C21_ = dcmplxMat(NZ_,NZ_);
C22_ = dcmplxMat(NZ_,NZ_);
C23_ = dcmplxMat(NZ_,NZ_);
C24_ = dcmplxMat(NZ_,NZ_);
C31_ = dcmplxMat(NZ_,NZ_);
C32_ = dcmplxMat(NZ_,NZ_);
C33_ = dcmplxMat(NZ_,NZ_);
C34_ = dcmplxMat(NZ_,NZ_);
C41_ = dcmplxMat(NZ_,NZ_);
C42_ = dcmplxMat(NZ_,NZ_);
C43_ = dcmplxMat(NZ_,NZ_);
C44_ = dcmplxMat(NZ_,NZ_);








//Assign automatically generated variables in equations 
// Submatrix variable definition (automatic)
C11_ = Ckl_in[i].block( 0, 0, NZ_, NZ_);
C12_ = Ckl_in[i].block( 0, NZ_, NZ_, NZ_);
C13_ = Ckl_in[i].block( 0, 2*NZ_, NZ_, NZ_);
C14_ = Ckl_in[i].block( 0, 3*NZ_, NZ_, NZ_);
C21_ = Ckl_in[i].block( NZ_, 0, NZ_, NZ_);
C22_ = Ckl_in[i].block( NZ_, NZ_, NZ_, NZ_);
C23_ = Ckl_in[i].block( NZ_, 2*NZ_, NZ_, NZ_);
C24_ = Ckl_in[i].block( NZ_, 3*NZ_, NZ_, NZ_);
C31_ = Ckl_in[i].block( 2*NZ_, 0, NZ_, NZ_);
C32_ = Ckl_in[i].block( 2*NZ_, NZ_, NZ_, NZ_);
C33_ = Ckl_in[i].block( 2*NZ_, 2*NZ_, NZ_, NZ_);
C34_ = Ckl_in[i].block( 2*NZ_, 3*NZ_, NZ_, NZ_);
C41_ = Ckl_in[i].block( 3*NZ_, 0, NZ_, NZ_);
C42_ = Ckl_in[i].block( 3*NZ_, NZ_, NZ_, NZ_);
C43_ = Ckl_in[i].block( 3*NZ_, 2*NZ_, NZ_, NZ_);
C44_ = Ckl_in[i].block( 3*NZ_, 3*NZ_, NZ_, NZ_);

// Scalar variable definition (automatic)
Tky = (kyctmp_)*fft2Dfac_;
Tm2TkxbTkyTq = ((-2.*(kxctmp_*kyctmp_))*q_)*fft2Dfac_;
Tmkxb = ((-1.*kxctmp_))*fft2Dfac_;

// Vector variable definition (automatic)
T2Tdz = (2.*kz_.matrix())*fft2Dfac_;
T2TdzTiLap2TkxbP2Tky = ((2.*(kyctmp_*pow(kxctmp_,2)))*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
T2TdzTiLap2Tky = ((2.*kyctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
T2TiLap2TkxbTkyP2 = ((2.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
TB0zTdz = (B0z_*kz_.matrix())*fft2Dfac_;
TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = ((-1.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix() + (ilap2tmp_*kz2_).matrix()*kxctmp_)*fft2Dfac_;
TdzTiLap2Tkxb = ((ilap2tmp_*kz_).matrix()*kxctmp_)*fft2Dfac_;
TdzTqPLUSTm2Tdz = (-2.*kz_.matrix() + kz_.matrix()*q_)*fft2Dfac_;
TiLap2Tky = (ilap2tmp_.matrix()*kyctmp_)*fft2Dfac_;
TmdzTiLap2Tkxb = ((-1.*kxctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
TmdzTq = ((-1.*q_)*kz_.matrix())*fft2Dfac_;
TmiLap2Tky = ((-1.*kyctmp_)*ilap2tmp_.matrix())*fft2Dfac_;







////////////////////////////////////////////////////////
///       AUTOMATICALLY GENERATED EQUATIONS         /////
///       see GenerateC++Equations.nb in MMA        /////

Ctmp_1_ = Tky*C31_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C31_ + T2TiLap2TkxbTkyP2.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
fft_.for_2DFull( Ctmp_1_ );
fft_.for_2DFull( Ctmp_2_ );
Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + (T2Tdz/fft2Dfac_).asDiagonal()*C21_ + (Tm2TkxbTkyTq/fft2Dfac_)*C11_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C31_;
dealias(Ctmp_2_);
Ckl_out[i].block( 0, 0, NZ_, NZ_) = Ctmp_2_;



Ctmp_1_ = Tky*C32_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C32_ + T2TiLap2TkxbTkyP2.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
fft_.for_2DFull( Ctmp_1_ );
fft_.for_2DFull( Ctmp_2_ );
Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + (T2Tdz/fft2Dfac_).asDiagonal()*C22_ + (Tm2TkxbTkyTq/fft2Dfac_)*C12_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C32_;
dealias(Ctmp_2_);
Ckl_out[i].block( 0, NZ_, NZ_, NZ_) = Ctmp_2_;



Ctmp_1_ = Tky*C33_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C33_ + T2TiLap2TkxbTkyP2.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
fft_.for_2DFull( Ctmp_1_ );
fft_.for_2DFull( Ctmp_2_ );
Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + (T2Tdz/fft2Dfac_).asDiagonal()*C23_ + (Tm2TkxbTkyTq/fft2Dfac_)*C13_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C33_;
dealias(Ctmp_2_);
Ckl_out[i].block( 0, 2*NZ_, NZ_, NZ_) = Ctmp_2_;



Ctmp_1_ = Tky*C34_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C34_ + T2TiLap2TkxbTkyP2.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
fft_.for_2DFull( Ctmp_1_ );
fft_.for_2DFull( Ctmp_2_ );
Ctmp_2_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_2_ + (T2Tdz/fft2Dfac_).asDiagonal()*C24_ + (Tm2TkxbTkyTq/fft2Dfac_)*C14_) + (TB0zTdz/fft2Dfac_).asDiagonal()*C34_;
dealias(Ctmp_2_);
Ckl_out[i].block( 0, 3*NZ_, NZ_, NZ_) = Ctmp_2_;



Ctmp_1_ = Tky*C41_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C31_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C41_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C11_;
dealias(Ctmp_3_);
Ckl_out[i].block( NZ_, 0, NZ_, NZ_) = Ctmp_3_;



Ctmp_1_ = Tky*C42_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C32_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C42_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C12_;
dealias(Ctmp_3_);
Ckl_out[i].block( NZ_, NZ_, NZ_, NZ_) = Ctmp_3_;



Ctmp_1_ = Tky*C43_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C33_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C43_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C13_;
dealias(Ctmp_3_);
Ckl_out[i].block( NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_3_;



Ctmp_1_ = Tky*C44_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzdzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C34_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C44_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C14_;
dealias(Ctmp_3_);
Ckl_out[i].block( NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_3_;



Ctmp_1_ = Tky*C11_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_1_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C11_;
dealias(Ctmp_1_);
Ckl_out[i].block( 2*NZ_, 0, NZ_, NZ_) = Ctmp_1_;



Ctmp_1_ = Tky*C12_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_1_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C12_;
dealias(Ctmp_1_);
Ckl_out[i].block( 2*NZ_, NZ_, NZ_, NZ_) = Ctmp_1_;



Ctmp_1_ = Tky*C13_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_1_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C13_;
dealias(Ctmp_1_);
Ckl_out[i].block( 2*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_1_;



Ctmp_1_ = Tky*C14_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_1_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C14_;
dealias(Ctmp_1_);
Ckl_out[i].block( 2*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_1_;



Ctmp_1_ = Tky*C21_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2Tky.asDiagonal()*C21_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C11_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C21_ + (TmdzTq/fft2Dfac_).asDiagonal()*C31_;
dealias(Ctmp_3_);
Ckl_out[i].block( 3*NZ_, 0, NZ_, NZ_) = Ctmp_3_;



Ctmp_1_ = Tky*C22_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2Tky.asDiagonal()*C22_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C12_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C22_ + (TmdzTq/fft2Dfac_).asDiagonal()*C32_;
dealias(Ctmp_3_);
Ckl_out[i].block( 3*NZ_, NZ_, NZ_, NZ_) = Ctmp_3_;



Ctmp_1_ = Tky*C23_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2Tky.asDiagonal()*C23_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C13_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C23_ + (TmdzTq/fft2Dfac_).asDiagonal()*C33_;
dealias(Ctmp_3_);
Ckl_out[i].block( 3*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_3_;



Ctmp_1_ = Tky*C24_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = By_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = T2TdzTiLap2Tky.asDiagonal()*C24_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C14_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = dzBy_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = dzdzBy_.asDiagonal()*Ctmp_3_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_3_ = Ctmp_1_ + (TB0zTdz/fft2Dfac_).asDiagonal()*C24_ + (TmdzTq/fft2Dfac_).asDiagonal()*C34_;
dealias(Ctmp_3_);
Ckl_out[i].block( 3*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_3_;
///     AUTOMATICALLY GENERATED EQUATIONS - END       ///
////////////////////////////////////////////////////////
