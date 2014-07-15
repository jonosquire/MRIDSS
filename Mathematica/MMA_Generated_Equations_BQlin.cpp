// Automatically generated temporary variables - class definition
dcmplxVecM T2Tdz, T2TdzTiLap2TkxbP2Tky, T2TdzTiLap2Tky, T2TiLap2TkxbTkyP2, TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2, TdzTiLap2Tkxb, TdzTqPLUSTm2Tdz, TiLap2Tky, TmdzTiLap2Tkxb, TmdzTq, TmiLap2Tky;

dcmplx Tky, Tm2TkxbTkyTq, Tmkxb;






// Automatically generated temporary variables - class constructor
T2Tdz = dcmplxVecM(NZ_);
T2TdzTiLap2TkxbP2Tky = dcmplxVecM(NZ_);
T2TdzTiLap2Tky = dcmplxVecM(NZ_);
T2TiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = dcmplxVecM(NZ_);
TdzTiLap2Tkxb = dcmplxVecM(NZ_);
TdzTqPLUSTm2Tdz = dcmplxVecM(NZ_);
TiLap2Tky = dcmplxVecM(NZ_);
TmdzTiLap2Tkxb = dcmplxVecM(NZ_);
TmdzTq = dcmplxVecM(NZ_);
TmiLap2Tky = dcmplxVecM(NZ_);








//Assign automatically generated variables in equations 
// Scalar variable definition (automatic)
Tky = kyctmp_;
Tm2TkxbTkyTq = (-2.*kxctmp_*kyctmp_)*q_;
Tmkxb = (-1.*kxctmp_);

// Vector variable definition (automatic)
T2Tdz = 2.*kz_.matrix();
T2TdzTiLap2TkxbP2Tky = (2.*kxctmp_*kxctmp_*kyctmp_)*(ilap2tmp_*kz_).matrix();
T2TdzTiLap2Tky = (2.*kyctmp_)*(ilap2tmp_*kz_).matrix();
T2TiLap2TkxbTkyP2 = (2.*kxctmp_*kyctmp_*kyctmp_)*ilap2tmp_.matrix();
TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2 = (-1.*kxctmp_*kyctmp_*kyctmp_)*ilap2tmp_.matrix() + (ilap2tmp_*kz_*kz_).matrix()*kxctmp_;
TdzTiLap2Tkxb = (ilap2tmp_*kz_).matrix()*kxctmp_;
TdzTqPLUSTm2Tdz = -2.*kz_.matrix() + kz_.matrix()*q_;
TiLap2Tky = ilap2tmp_.matrix()*kyctmp_;
TmdzTiLap2Tkxb = (-1.*kxctmp_)*(ilap2tmp_*kz_).matrix();
TmdzTq = (-1.*q_)*kz_.matrix();
TmiLap2Tky = (-1.*kyctmp_)*ilap2tmp_.matrix();







////////////////////////////////////////////////////////
///       AUTOMATICALLY GENERATED EQUATIONS         /////
///       see GenerateC++Equations.nb in MMA        /////

Ckl_out[i].block( 0, 0, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C31_ )))+ilapFtmp_.matrix().asDiagonal()*(fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2TkxbP2Tky.asDiagonal()*C31_ )+fft_.back_2DFull(T2TiLap2TkxbTkyP2.asDiagonal()*C41_ ))))+((T2Tdz.asDiagonal()*C21_ )+(Tm2TkxbTkyTq*C11_ )))  );


Ckl_out[i].block( 0, NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C32_ )))+ilapFtmp_.matrix().asDiagonal()*(fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2TkxbP2Tky.asDiagonal()*C32_ )+fft_.back_2DFull(T2TiLap2TkxbTkyP2.asDiagonal()*C42_ ))))+((T2Tdz.asDiagonal()*C22_ )+(Tm2TkxbTkyTq*C12_ )))  );


Ckl_out[i].block( 0, 2*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C33_ )))+ilapFtmp_.matrix().asDiagonal()*(fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2TkxbP2Tky.asDiagonal()*C33_ )+fft_.back_2DFull(T2TiLap2TkxbTkyP2.asDiagonal()*C43_ ))))+((T2Tdz.asDiagonal()*C23_ )+(Tm2TkxbTkyTq*C13_ )))  );


Ckl_out[i].block( 0, 3*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C34_ )))+ilapFtmp_.matrix().asDiagonal()*(fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2TkxbP2Tky.asDiagonal()*C34_ )+fft_.back_2DFull(T2TiLap2TkxbTkyP2.asDiagonal()*C44_ ))))+((T2Tdz.asDiagonal()*C24_ )+(Tm2TkxbTkyTq*C14_ )))  );


Ckl_out[i].block( NZ_, 0, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C41_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*(fft_.back_2DFull(Tmkxb*C31_ )))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TmdzTiLap2Tkxb.asDiagonal()*C31_ )+fft_.back_2DFull(TmiLap2Tky.asDiagonal()*C41_ ))))+(TdzTqPLUSTm2Tdz.asDiagonal()*C11_ )  );


Ckl_out[i].block( NZ_, NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C42_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*(fft_.back_2DFull(Tmkxb*C32_ )))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TmdzTiLap2Tkxb.asDiagonal()*C32_ )+fft_.back_2DFull(TmiLap2Tky.asDiagonal()*C42_ ))))+(TdzTqPLUSTm2Tdz.asDiagonal()*C12_ )  );


Ckl_out[i].block( NZ_, 2*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C43_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*(fft_.back_2DFull(Tmkxb*C33_ )))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TmdzTiLap2Tkxb.asDiagonal()*C33_ )+fft_.back_2DFull(TmiLap2Tky.asDiagonal()*C43_ ))))+(TdzTqPLUSTm2Tdz.asDiagonal()*C13_ )  );


Ckl_out[i].block( NZ_, 3*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C44_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*(fft_.back_2DFull(Tmkxb*C34_ )))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TmdzTiLap2Tkxb.asDiagonal()*C34_ )+fft_.back_2DFull(TmiLap2Tky.asDiagonal()*C44_ ))))+(TdzTqPLUSTm2Tdz.asDiagonal()*C14_ )  );


Ckl_out[i].block( 2*NZ_, 0, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C11_ )))  );


Ckl_out[i].block( 2*NZ_, NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C12_ )))  );


Ckl_out[i].block( 2*NZ_, 2*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C13_ )))  );


Ckl_out[i].block( 2*NZ_, 3*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C14_ )))  );


Ckl_out[i].block( 3*NZ_, 0, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C21_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2Tky.asDiagonal()*C21_ )+fft_.back_2DFull(TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C11_ ))))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TdzTiLap2Tkxb.asDiagonal()*C11_ )+fft_.back_2DFull(TiLap2Tky.asDiagonal()*C21_ ))))+(TmdzTq.asDiagonal()*C31_ )  );


Ckl_out[i].block( 3*NZ_, NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C22_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2Tky.asDiagonal()*C22_ )+fft_.back_2DFull(TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C12_ ))))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TdzTiLap2Tkxb.asDiagonal()*C12_ )+fft_.back_2DFull(TiLap2Tky.asDiagonal()*C22_ ))))+(TmdzTq.asDiagonal()*C32_ )  );


Ckl_out[i].block( 3*NZ_, 2*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C23_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2Tky.asDiagonal()*C23_ )+fft_.back_2DFull(TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C13_ ))))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TdzTiLap2Tkxb.asDiagonal()*C13_ )+fft_.back_2DFull(TiLap2Tky.asDiagonal()*C23_ ))))+(TmdzTq.asDiagonal()*C33_ )  );


Ckl_out[i].block( 3*NZ_, 3*NZ_, NZ_, NZ_)  =  dealias(  fft_.for_2DFull(By_.asDiagonal()*(fft_.back_2DFull(Tky*C24_ )))+fft_.for_2DFull(dzBy_.asDiagonal()*((fft_.back_2DFull(T2TdzTiLap2Tky.asDiagonal()*C24_ )+fft_.back_2DFull(TdzP2TiLap2TkxbPLUSTmiLap2TkxbTkyP2.asDiagonal()*C14_ ))))+fft_.for_2DFull(dzdzBy_.asDiagonal()*((fft_.back_2DFull(TdzTiLap2Tkxb.asDiagonal()*C14_ )+fft_.back_2DFull(TiLap2Tky.asDiagonal()*C24_ ))))+(TmdzTq.asDiagonal()*C34_ )  );


///     AUTOMATICALLY GENERATED EQUATIONS - END       ///
////////////////////////////////////////////////////////
