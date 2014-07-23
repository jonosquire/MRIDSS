// Automatically generated temporary variables - class definition
dcmplxVecM T2Tdz, T2TdzP2TiLap2Tkxb, T2TdzTiLap2TkxbP2Tky, T2TdzTiLap2Tky, T2TiLap2TkxbTkyP2, TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2, TdzTiLap2Tkxb, TdzTiLap2TkxbP3PLUSTmdzTkxb, TdzTkxbPLUSTmdzTiLap2TkxbP3, TdzTqPLUSTm2Tdz, TiLap2TkxbP2Tky, TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb, TiLap2Tky, Tm2TdzP2TiLap2Tkxb, Tm2TdzTiLap2TkxbP2Tky, Tm2TdzTiLap2Tky, Tm2TiLap2TkxbTkyP2, TmdzTiLap2Tkxb, TmdzTq, TmiLap2TkxbP2Tky, TmiLap2Tky;



dcmplx Tkxb, Tky, Tm2TkxbTkyTq, Tmkxb, Tmky;

dcmplxMat Ctmp_1_, Ctmp_2_, Ctmp_3_, Ctmp_4_, Ctmp_5_, Ctmp_6_, Ctmp_7_, Ctmp_8_, Ctmp_9_, Ctmp_10_, Ctmp_11_, Ctmp_12_;

dcmplxMat C11_, C12_, C13_, C14_, C21_, C22_, C23_, C24_, C31_, C32_, C33_, C34_, C41_, C42_, C43_, C44_;






// Automatically generated temporary variables - class constructor
T2Tdz = dcmplxVecM(NZ_);
T2TdzP2TiLap2Tkxb = dcmplxVecM(NZ_);
T2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
T2TdzTiLap2Tky = dcmplxVecM(NZ_);
T2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TdzTiLap2Tkxb = dcmplxVecM(NZ_);
TdzTiLap2TkxbP3PLUSTmdzTkxb = dcmplxVecM(NZ_);
TdzTkxbPLUSTmdzTiLap2TkxbP3 = dcmplxVecM(NZ_);
TdzTqPLUSTm2Tdz = dcmplxVecM(NZ_);
TiLap2TkxbP2Tky = dcmplxVecM(NZ_);
TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb = dcmplxVecM(NZ_);
TiLap2Tky = dcmplxVecM(NZ_);
Tm2TdzP2TiLap2Tkxb = dcmplxVecM(NZ_);
Tm2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
Tm2TdzTiLap2Tky = dcmplxVecM(NZ_);
Tm2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TmdzTiLap2Tkxb = dcmplxVecM(NZ_);
TmdzTq = dcmplxVecM(NZ_);
TmiLap2TkxbP2Tky = dcmplxVecM(NZ_);
TmiLap2Tky = dcmplxVecM(NZ_);

Ctmp_1_ = dcmplxMat(NZ_,NZ_);
Ctmp_2_ = dcmplxMat(NZ_,NZ_);
Ctmp_3_ = dcmplxMat(NZ_,NZ_);
Ctmp_4_ = dcmplxMat(NZ_,NZ_);
Ctmp_5_ = dcmplxMat(NZ_,NZ_);
Ctmp_6_ = dcmplxMat(NZ_,NZ_);
Ctmp_7_ = dcmplxMat(NZ_,NZ_);
Ctmp_8_ = dcmplxMat(NZ_,NZ_);
Ctmp_9_ = dcmplxMat(NZ_,NZ_);
Ctmp_10_ = dcmplxMat(NZ_,NZ_);
Ctmp_11_ = dcmplxMat(NZ_,NZ_);
Ctmp_12_ = dcmplxMat(NZ_,NZ_);

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
Tkxb = (kxctmp_)*fft2Dfac_;
Tky = (kyctmp_)*fft2Dfac_;
Tm2TkxbTkyTq = ((-2.*(kxctmp_*kyctmp_))*q_)*fft2Dfac_;
Tmkxb = ((-1.*kxctmp_))*fft2Dfac_;
Tmky = ((-1.*kyctmp_))*fft2Dfac_;

// Vector variable definition (automatic)
T2Tdz = (2.*kz_.matrix())*fft2Dfac_;
T2TdzP2TiLap2Tkxb = ((2.*kxctmp_)*(ilap2tmp_*kz2_).matrix())*fft2Dfac_;
T2TdzTiLap2TkxbP2Tky = ((2.*(kyctmp_*pow(kxctmp_,2)))*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
T2TdzTiLap2Tky = ((2.*kyctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
T2TiLap2TkxbTkyP2 = ((2.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = ((-1.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix() + (ilap2tmp_*kz2_).matrix()*kxctmp_)*fft2Dfac_;
TdzTiLap2Tkxb = ((ilap2tmp_*kz_).matrix()*kxctmp_)*fft2Dfac_;
TdzTiLap2TkxbP3PLUSTmdzTkxb = ((-1.*kxctmp_)*kz_.matrix() + (ilap2tmp_*kz_).matrix()*pow(kxctmp_,3))*fft2Dfac_;
TdzTkxbPLUSTmdzTiLap2TkxbP3 = ((-1.*pow(kxctmp_,3))*(ilap2tmp_*kz_).matrix() + kxctmp_*kz_.matrix())*fft2Dfac_;
TdzTqPLUSTm2Tdz = (-2.*kz_.matrix() + kz_.matrix()*q_)*fft2Dfac_;
TiLap2TkxbP2Tky = (ilap2tmp_.matrix()*(kyctmp_*pow(kxctmp_,2)))*fft2Dfac_;
TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb = ((-1.*kxctmp_)*(ilap2tmp_*kz2_).matrix() + ilap2tmp_.matrix()*(kxctmp_*pow(kyctmp_,2)))*fft2Dfac_;
TiLap2Tky = (ilap2tmp_.matrix()*kyctmp_)*fft2Dfac_;
Tm2TdzP2TiLap2Tkxb = ((-2.*kxctmp_)*(ilap2tmp_*kz2_).matrix())*fft2Dfac_;
Tm2TdzTiLap2TkxbP2Tky = ((-2.*(kyctmp_*pow(kxctmp_,2)))*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
Tm2TdzTiLap2Tky = ((-2.*kyctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
Tm2TiLap2TkxbTkyP2 = ((-2.*(kxctmp_*pow(kyctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
TmdzTiLap2Tkxb = ((-1.*kxctmp_)*(ilap2tmp_*kz_).matrix())*fft2Dfac_;
TmdzTq = ((-1.*q_)*kz_.matrix())*fft2Dfac_;
TmiLap2TkxbP2Tky = ((-1.*(kyctmp_*pow(kxctmp_,2)))*ilap2tmp_.matrix())*fft2Dfac_;
TmiLap2Tky = ((-1.*kyctmp_)*ilap2tmp_.matrix())*fft2Dfac_;







////////////////////////////////////////////////////////
///       AUTOMATICALLY GENERATED EQUATIONS         /////
///       see GenerateC++Equations.nb in MMA        /////

Ctmp_1_ = Tkxb*C31_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C31_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C11_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C11_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C11_ + T2TdzTiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C31_ + T2TiLap2TkxbTkyP2.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C31_ + Tm2TdzTiLap2Tky.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C11_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C31_ + TiLap2TkxbP2Tky.asDiagonal()*C41_ + Tmky*C41_;
fft_.back_2DFull( Ctmp_11_ );
Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C11_ + Tky*C21_ + TmiLap2TkxbP2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_12_ );
Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_5_ );
Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C21_ + (Tm2TkxbTkyTq/fft2Dfac_)*C11_);
dealias(Ctmp_12_);
Ckl_out[i].block( 0, 0, NZ_, NZ_) = Ctmp_12_;



Ctmp_1_ = Tkxb*C32_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C32_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C12_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C12_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C12_ + T2TdzTiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C32_ + T2TiLap2TkxbTkyP2.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C32_ + Tm2TdzTiLap2Tky.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C12_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C32_ + TiLap2TkxbP2Tky.asDiagonal()*C42_ + Tmky*C42_;
fft_.back_2DFull( Ctmp_11_ );
Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C12_ + Tky*C22_ + TmiLap2TkxbP2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_12_ );
Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_5_ );
Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C22_ + (Tm2TkxbTkyTq/fft2Dfac_)*C12_);
dealias(Ctmp_12_);
Ckl_out[i].block( 0, NZ_, NZ_, NZ_) = Ctmp_12_;



Ctmp_1_ = Tkxb*C33_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C33_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C13_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C13_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C13_ + T2TdzTiLap2Tky.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C33_ + T2TiLap2TkxbTkyP2.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C33_ + Tm2TdzTiLap2Tky.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C13_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C33_ + TiLap2TkxbP2Tky.asDiagonal()*C43_ + Tmky*C43_;
fft_.back_2DFull( Ctmp_11_ );
Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C13_ + Tky*C23_ + TmiLap2TkxbP2Tky.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_12_ );
Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_5_ );
Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C23_ + (Tm2TkxbTkyTq/fft2Dfac_)*C13_);
dealias(Ctmp_12_);
Ckl_out[i].block( 0, 2*NZ_, NZ_, NZ_) = Ctmp_12_;



Ctmp_1_ = Tkxb*C34_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C34_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C14_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C14_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzP2TiLap2Tkxb.asDiagonal()*C14_ + T2TdzTiLap2Tky.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = T2TdzTiLap2TkxbP2Tky.asDiagonal()*C34_ + T2TiLap2TkxbTkyP2.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzdzdzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = Tm2TdzP2TiLap2Tkxb.asDiagonal()*C34_ + Tm2TdzTiLap2Tky.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tm2TdzTiLap2TkxbP2Tky.asDiagonal()*C14_ + Tm2TiLap2TkxbTkyP2.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzdzBx_.asDiagonal()*Ctmp_10_;
Ctmp_11_ = TdzTiLap2TkxbP3PLUSTmdzTkxb.asDiagonal()*C34_ + TiLap2TkxbP2Tky.asDiagonal()*C44_ + Tmky*C44_;
fft_.back_2DFull( Ctmp_11_ );
Ctmp_11_ = dzBx_.asDiagonal()*Ctmp_11_;
Ctmp_12_ = TdzTkxbPLUSTmdzTiLap2TkxbP3.asDiagonal()*C14_ + Tky*C24_ + TmiLap2TkxbP2Tky.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_12_ );
Ctmp_12_ = dzUx_.asDiagonal()*Ctmp_12_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_5_ = Ctmp_10_ + Ctmp_11_ + Ctmp_12_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_5_ );
Ctmp_12_ = Ctmp_1_ + ilapFtmp_.matrix().asDiagonal()*(Ctmp_5_ + (T2Tdz/fft2Dfac_).asDiagonal()*C24_ + (Tm2TkxbTkyTq/fft2Dfac_)*C14_);
dealias(Ctmp_12_);
Ckl_out[i].block( 0, 3*NZ_, NZ_, NZ_) = Ctmp_12_;



Ctmp_1_ = Tkxb*C41_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C41_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C21_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C21_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C41_ + TmiLap2TkxbP2Tky.asDiagonal()*C31_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C11_ + TmdzTiLap2Tkxb.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tkxb*C11_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = Tmkxb*C31_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C11_;
dealias(Ctmp_10_);
Ckl_out[i].block( NZ_, 0, NZ_, NZ_) = Ctmp_10_;



Ctmp_1_ = Tkxb*C42_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C42_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C22_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C22_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C42_ + TmiLap2TkxbP2Tky.asDiagonal()*C32_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C12_ + TmdzTiLap2Tkxb.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tkxb*C12_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = Tmkxb*C32_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C12_;
dealias(Ctmp_10_);
Ckl_out[i].block( NZ_, NZ_, NZ_, NZ_) = Ctmp_10_;



Ctmp_1_ = Tkxb*C43_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C43_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C23_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C23_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C43_ + TmiLap2TkxbP2Tky.asDiagonal()*C33_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C13_ + TmdzTiLap2Tkxb.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tkxb*C13_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = Tmkxb*C33_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C13_;
dealias(Ctmp_10_);
Ckl_out[i].block( NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_10_;



Ctmp_1_ = Tkxb*C44_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C44_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C24_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C24_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzdzUy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C44_ + TmiLap2TkxbP2Tky.asDiagonal()*C34_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzBx_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TiLap2TkxbP2Tky.asDiagonal()*C14_ + TmdzTiLap2Tkxb.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzUx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzdzBy_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = Tkxb*C14_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = Tmkxb*C34_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzBy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TdzTqPLUSTm2Tdz/fft2Dfac_).asDiagonal()*C14_;
dealias(Ctmp_10_);
Ckl_out[i].block( NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_10_;



Ctmp_1_ = Tkxb*C11_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C11_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C31_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C31_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_6_ = Ctmp_1_;
dealias(Ctmp_6_);
Ckl_out[i].block( 2*NZ_, 0, NZ_, NZ_) = Ctmp_6_;



Ctmp_1_ = Tkxb*C12_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C12_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C32_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C32_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_6_ = Ctmp_1_;
dealias(Ctmp_6_);
Ckl_out[i].block( 2*NZ_, NZ_, NZ_, NZ_) = Ctmp_6_;



Ctmp_1_ = Tkxb*C13_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C13_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C33_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C33_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_6_ = Ctmp_1_;
dealias(Ctmp_6_);
Ckl_out[i].block( 2*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_6_;



Ctmp_1_ = Tkxb*C14_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C14_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C34_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C34_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBx_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzUx_.asDiagonal()*Ctmp_6_;
Ctmp_1_ = Ctmp_1_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_6_ = Ctmp_1_;
dealias(Ctmp_6_);
Ckl_out[i].block( 2*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_6_;



Ctmp_1_ = Tkxb*C21_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C21_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C41_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C41_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C21_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C11_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C11_ + TiLap2Tky.asDiagonal()*C21_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C21_ + TmiLap2TkxbP2Tky.asDiagonal()*C11_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C31_ + TmdzTiLap2Tkxb.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C31_ + Tm2TdzTiLap2Tky.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C31_ + TmiLap2Tky.asDiagonal()*C41_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TmdzTq/fft2Dfac_).asDiagonal()*C31_;
dealias(Ctmp_10_);
Ckl_out[i].block( 3*NZ_, 0, NZ_, NZ_) = Ctmp_10_;



Ctmp_1_ = Tkxb*C22_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C22_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C42_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C42_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C22_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C12_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C12_ + TiLap2Tky.asDiagonal()*C22_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C22_ + TmiLap2TkxbP2Tky.asDiagonal()*C12_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C32_ + TmdzTiLap2Tkxb.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C32_ + Tm2TdzTiLap2Tky.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C32_ + TmiLap2Tky.asDiagonal()*C42_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TmdzTq/fft2Dfac_).asDiagonal()*C32_;
dealias(Ctmp_10_);
Ckl_out[i].block( 3*NZ_, NZ_, NZ_, NZ_) = Ctmp_10_;



Ctmp_1_ = Tkxb*C23_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C23_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C43_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C43_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C23_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C13_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C13_ + TiLap2Tky.asDiagonal()*C23_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C23_ + TmiLap2TkxbP2Tky.asDiagonal()*C13_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C33_ + TmdzTiLap2Tkxb.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C33_ + Tm2TdzTiLap2Tky.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C33_ + TmiLap2Tky.asDiagonal()*C43_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TmdzTq/fft2Dfac_).asDiagonal()*C33_;
dealias(Ctmp_10_);
Ckl_out[i].block( 3*NZ_, 2*NZ_, NZ_, NZ_) = Ctmp_10_;



Ctmp_1_ = Tkxb*C24_;
fft_.back_2DFull( Ctmp_1_ );
Ctmp_1_ = Bx_.asDiagonal()*Ctmp_1_;
Ctmp_2_ = Tky*C24_;
fft_.back_2DFull( Ctmp_2_ );
Ctmp_2_ = By_.asDiagonal()*Ctmp_2_;
Ctmp_3_ = Tmkxb*C44_;
fft_.back_2DFull( Ctmp_3_ );
Ctmp_3_ = Ux_.asDiagonal()*Ctmp_3_;
Ctmp_4_ = Tmky*C44_;
fft_.back_2DFull( Ctmp_4_ );
Ctmp_4_ = Uy_.asDiagonal()*Ctmp_4_;
Ctmp_5_ = T2TdzTiLap2Tky.asDiagonal()*C24_ + TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C14_;
fft_.back_2DFull( Ctmp_5_ );
Ctmp_5_ = dzBy_.asDiagonal()*Ctmp_5_;
Ctmp_6_ = TdzTiLap2Tkxb.asDiagonal()*C14_ + TiLap2Tky.asDiagonal()*C24_;
fft_.back_2DFull( Ctmp_6_ );
Ctmp_6_ = dzdzBy_.asDiagonal()*Ctmp_6_;
Ctmp_7_ = TdzTiLap2Tkxb.asDiagonal()*C24_ + TmiLap2TkxbP2Tky.asDiagonal()*C14_;
fft_.back_2DFull( Ctmp_7_ );
Ctmp_7_ = dzBx_.asDiagonal()*Ctmp_7_;
Ctmp_8_ = TiLap2TkxbP2Tky.asDiagonal()*C34_ + TmdzTiLap2Tkxb.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_8_ );
Ctmp_8_ = dzUx_.asDiagonal()*Ctmp_8_;
Ctmp_9_ = TiLap2TkxbTkyP2PLUSTmdzP2TiLap2Tkxb.asDiagonal()*C34_ + Tm2TdzTiLap2Tky.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_9_ );
Ctmp_9_ = dzUy_.asDiagonal()*Ctmp_9_;
Ctmp_10_ = TmdzTiLap2Tkxb.asDiagonal()*C34_ + TmiLap2Tky.asDiagonal()*C44_;
fft_.back_2DFull( Ctmp_10_ );
Ctmp_10_ = dzdzUy_.asDiagonal()*Ctmp_10_;
Ctmp_1_ = Ctmp_1_ + Ctmp_10_ + Ctmp_2_ + Ctmp_3_ + Ctmp_4_ + Ctmp_5_ + Ctmp_6_ + Ctmp_7_ + Ctmp_8_ + Ctmp_9_;
fft_.for_2DFull( Ctmp_1_ );
Ctmp_10_ = Ctmp_1_ + (TmdzTq/fft2Dfac_).asDiagonal()*C34_;
dealias(Ctmp_10_);
Ckl_out[i].block( 3*NZ_, 3*NZ_, NZ_, NZ_) = Ctmp_10_;
///     AUTOMATICALLY GENERATED EQUATIONS - END       ///
////////////////////////////////////////////////////////
