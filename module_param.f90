module variables

  use decomp_2d, only : mytype

  ! Boundary conditions : ncl = 2 --> Dirichlet
  ! Boundary conditions : ncl = 1 --> Free-slip
  ! Boundary conditions : ncl = 0 --> Periodic
  ! l: power of 2,3,4,5 and 6
  ! if ncl = 1 or 2, --> n  = 2l+ 1
  !                  --> nm = n - 1
  !                  --> m  = n + 1
  ! If ncl = 0,      --> n  = 2*l
  !                  --> nm = n
  !                  --> m  = n + 2
  !nstat = size arrays for statistic collection
  !2-->every 2 mesh nodes
  !4-->every 4 mesh nodes
  !nvisu = size for visualization collection
  !nprobe =  size for probe collection (energy spectra)

!Possible n points: 3 5 7 9 11 13 17 19 21 25 31 33 37 41 49 51 55 61 65 73 81 91 97 101 109 121 129 145 151 161 163 181 193 201 217 241 251 257 271 289 301 321 325 361 385 401 433 451 481 487 501 513 541 577 601 641 649 721 751 769 801 811 865 901 961 973 1001 1025 1081 1153 1201 1251 1281 1297 1351 1441 1459 1501 1537 1601 1621 1729 1801 1921 1945 2001 2049 2161 2251 2305 2401 2431 2501 2561 2593 2701 2881 2917 3001 3073 3201 3241 3457 3601 3751 3841 3889 4001 4051 4097 4321 4375 4501 4609 4801 4861 5001 5121 5185 5401 5761 5833 6001 6145 6251 6401 6481 6751 6913 7201 7291 7501 7681 7777 8001 8101 8193 8641 8749 9001 9217 9601 9721 enough

  integer :: nx,ny,nz,nphi,p_row,p_col,nxm,nym,nzm
  integer :: nstat=1,nvisu=1,nprobe=1,nlength=1

  real(mytype),allocatable,dimension(:) :: nsc,uset,cp,ri,group
  real(mytype) :: fpi2

#ifndef DOUBLE_PREC
  integer,parameter :: prec = 4
#else
#ifdef SAVE_SINGLE
  integer,parameter :: prec = 4
#else
  integer,parameter :: prec = 8  
#endif
#endif

  !module filter
  real(mytype),allocatable,dimension(:) :: fifx,ficx,fibx,fiffx,fibbx,fiz1x,fiz2x
  real(mytype),allocatable,dimension(:,:) :: filax,filaxp
  real(mytype),allocatable,dimension(:) :: fifxp,ficxp,fibxp,fiffxp,fibbxp
  real(mytype),allocatable,dimension(:) :: fify,ficy,fiby,fiffy,fibby,fiz1y,fiz2y
  real(mytype),allocatable,dimension(:,:) :: filay,filayp
  real(mytype),allocatable,dimension(:) :: fifyp,ficyp,fibyp,fiffyp,fibbyp
  real(mytype),allocatable,dimension(:) :: fifz,ficz,fibz,fiffz,fibbz,fiz1z,fiz2z
  real(mytype),allocatable, dimension(:,:) :: filaz,filazp
  real(mytype),allocatable,dimension(:) :: fifzp,ficzp,fibzp,fiffzp,fibbzp
  integer,allocatable,dimension(:) :: idata
  real(mytype),allocatable,dimension(:) :: Cs

  !module derivative
  real(mytype),allocatable,dimension(:) :: ffx,sfx,fsx,fwx,ssx,swx
  real(mytype),allocatable,dimension(:) :: ffxp,sfxp,fsxp,fwxp,ssxp,swxp
  real(mytype),allocatable,dimension(:) :: ffy,sfy,fsy,fwy,ssy,swy
  real(mytype),allocatable,dimension(:) :: ffyp,sfyp,fsyp,fwyp,ssyp,swyp
  real(mytype),allocatable,dimension(:) :: ffz,sfz,fsz,fwz,ssz,swz
  real(mytype),allocatable,dimension(:) :: ffzp,sfzp,fszp,fwzp,sszp,swzp

  real(mytype),allocatable,dimension(:) :: ffxS,sfxS,fsxS,fwxS,ssxS,swxS
  real(mytype),allocatable,dimension(:) :: ffxpS,sfxpS,fsxpS,fwxpS,ssxpS,swxpS
  real(mytype),allocatable,dimension(:) :: ffyS,sfyS,fsyS,fwyS,ssyS,swyS
  real(mytype),allocatable,dimension(:) :: ffypS,sfypS,fsypS,fwypS,ssypS,swypS
  real(mytype),allocatable,dimension(:) :: ffzS,sfzS,fszS,fwzS,sszS,swzS
  real(mytype),allocatable,dimension(:) :: ffzpS,sfzpS,fszpS,fwzpS,sszpS,swzpS

  real(mytype), save, allocatable, dimension(:,:) :: sx,vx
  real(mytype), save, allocatable, dimension(:,:) :: sy,vy
  real(mytype), save, allocatable, dimension(:,:) :: sz,vz

  ABSTRACT INTERFACE
     SUBROUTINE DERIVATIVE_X(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r 
       real(mytype), dimension(ny,nz):: s
       real(mytype), dimension(nx):: ff,fs,fw 
     END SUBROUTINE DERIVATIVE_X
     SUBROUTINE DERIVATIVE_Y(t,u,r,s,ff,fs,fw,pp,nx,ny,nz,npaire)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r 
       real(mytype), dimension(nx,nz):: s
       real(mytype), dimension(ny):: ff,fs,fw,pp
     END SUBROUTINE DERIVATIVE_Y
     SUBROUTINE DERIVATIVE_YY(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r 
       real(mytype), dimension(nx,nz):: s
       real(mytype), dimension(ny):: ff,fs,fw
     END SUBROUTINE DERIVATIVE_YY
     SUBROUTINE DERIVATIVE_Z(t,u,r,s,ff,fs,fw,nx,ny,nz,npaire)
       use decomp_2d, only : mytype
       integer :: nx,ny,nz,npaire
       real(mytype), dimension(nx,ny,nz) :: t,u,r 
       real(mytype), dimension(nx,ny):: s
       real(mytype), dimension(nz):: ff,fs,fw 
     END SUBROUTINE DERIVATIVE_Z
  END INTERFACE

  PROCEDURE (DERIVATIVE_X) derx_00,derx_11,derx_12,derx_21,derx_22,&
       derxx_dns_00,derxx_dns_11,derxx_dns_12,derxx_dns_21,derxx_dns_22,&
       derxx_iles_00,derxx_iles_11,derxx_iles_12,derxx_iles_21,derxx_iles_22
  PROCEDURE (DERIVATIVE_X), POINTER :: derx,derxx,derxS,derxxS
  PROCEDURE (DERIVATIVE_Y) dery_00,dery_11,dery_12,dery_21,dery_22
  PROCEDURE (DERIVATIVE_Y), POINTER :: dery,deryS
  PROCEDURE (DERIVATIVE_YY) &
       deryy_dns_00,deryy_dns_11,deryy_dns_12,deryy_dns_21,deryy_dns_22,&
       deryy_iles_00,deryy_iles_11,deryy_iles_12,deryy_iles_21,deryy_iles_22
  PROCEDURE (DERIVATIVE_YY), POINTER :: deryy,deryyS
  PROCEDURE (DERIVATIVE_Z) derz_00,derz_11,derz_12,derz_21,derz_22,&
       derzz_dns_00,derzz_dns_11,derzz_dns_12,derzz_dns_21,derzz_dns_22,&
       derzz_iles_00,derzz_iles_11,derzz_iles_12,derzz_iles_21,derzz_iles_22
  PROCEDURE (DERIVATIVE_Z), POINTER :: derz,derzz,derzS,derzzS

  !O6SVV
  real(mytype),allocatable,dimension(:) :: newsm,newtm,newsmt,newtmt
  real(mytype),allocatable,dimension(:) :: newrm,ttm,newrmt,ttmt

  !module pressure
  real(mytype), save, allocatable, dimension(:,:) :: dpdyx1,dpdyxn,dpdzx1,dpdzxn
  real(mytype), save, allocatable, dimension(:,:) :: dpdxy1,dpdxyn,dpdzy1,dpdzyn
  real(mytype), save, allocatable, dimension(:,:) :: dpdxz1,dpdxzn,dpdyz1,dpdyzn

  !module inflow
  real(mytype), save, allocatable, dimension(:,:) :: bxx1,bxy1,bxz1,bxxn,bxyn,bxzn,bxo,byo,bzo
  real(mytype), save, allocatable, dimension(:,:) :: byx1,byy1,byz1,byxn,byyn,byzn
  real(mytype), save, allocatable, dimension(:,:) :: bzx1,bzy1,bzz1,bzxn,bzyn,bzzn

  !module derpres
  real(mytype),allocatable,dimension(:) :: cfx6,ccx6,cbx6,cfxp6,ciwxp6,csxp6,&
       cwxp6,csx6,cwx6,cifx6,cicx6,cisx6
  real(mytype),allocatable,dimension(:) :: cibx6,cifxp6,cisxp6,ciwx6
  real(mytype),allocatable,dimension(:) :: cfi6,cci6,cbi6,cfip6,csip6,cwip6,csi6,&
       cwi6,cifi6,cici6,cibi6,cifip6
  real(mytype),allocatable,dimension(:) :: cisip6,ciwip6,cisi6,ciwi6
  real(mytype),allocatable,dimension(:) :: cfy6,ccy6,cby6,cfyp6,csyp6,cwyp6,csy6
  real(mytype),allocatable,dimension(:) :: cwy6,cify6,cicy6,ciby6,cifyp6,cisyp6,&
       ciwyp6,cisy6,ciwy6
  real(mytype),allocatable,dimension(:) :: cfi6y,cci6y,cbi6y,cfip6y,csip6y,cwip6y,&
       csi6y,cwi6y,cifi6y,cici6y
  real(mytype),allocatable,dimension(:) :: cibi6y,cifip6y,cisip6y,ciwip6y,cisi6y,ciwi6y
  real(mytype),allocatable,dimension(:) :: cfz6,ccz6,cbz6,cfzp6,cszp6,cwzp6,csz6
  real(mytype),allocatable,dimension(:) :: cwz6,cifz6,cicz6,cibz6,cifzp6,ciszp6,&
       ciwzp6,cisz6,ciwz6
  real(mytype),allocatable,dimension(:) :: cfi6z,cci6z,cbi6z,cfip6z,csip6z,cwip6z,&
       csi6z,cwi6z,cifi6z,cici6z
  real(mytype),allocatable,dimension(:) :: cibi6z,cifip6z,cisip6z,ciwip6z,cisi6z,ciwi6z

  !module waves
  complex(mytype),allocatable,dimension(:) :: zkz,zk2,ezs
  complex(mytype),allocatable,dimension(:) :: yky,yk2,eys
  complex(mytype),allocatable,dimension(:) :: xkx,xk2,exs

  !module mesh
  real(mytype),allocatable,dimension(:) :: ppy,pp2y,pp4y
  real(mytype),allocatable,dimension(:) :: ppyi,pp2yi,pp4yi
  real(mytype),allocatable,dimension(:) :: yp,ypi,del
  real(mytype),allocatable,dimension(:) :: yeta,yetai
  real(mytype) :: alpha,beta

end module variables

module param

  use decomp_2d, only : mytype

  integer :: nclx1,nclxn,ncly1,nclyn,nclz1,nclzn
  integer :: nclxS1,nclxSn,nclyS1,nclySn,nclzS1,nclzSn

  !logical variable for boundary condition that is true in periodic case
  !and false otherwise
  logical :: nclx,ncly,nclz

  integer :: cont_phi,itr,itime
  integer :: ifft,ivirt,istret,iforc_entree,iturb
  integer :: itype,iin,nscheme,ifirst,ilast,iles
  integer :: isave,ilit,idebmod,imodulo,imodulo2,idemarre,icommence,irecord
  integer :: iscalar,nxboite,istat,iread,iadvance_time,irotation
  integer :: ilag,npif,izap
  real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2,t,xxk1,xxk2
  real(mytype) :: dt,re,xnu,noise,noise1,u1,u2,angle,anglex,angley
  real(mytype) :: wrotation

  !LES
  integer :: jLES
  real(mytype) :: smagcst,walecst,FSGS

  character :: filesauve*80, filenoise*80, &
  nchamp*80,filepath*80, fileturb*80, filevisu*80, datapath*80
  real(mytype), dimension(5) :: adt,bdt,cdt,gdt

  !VISU
  integer :: save_w,save_w1,save_w2,save_w3,save_qc,save_pc
  integer :: save_ux,save_uy,save_uz,save_phi,save_pre
  integer :: save_uxm,save_uym,save_uzm,save_phim,save_prem
  integer :: save_ibm,save_dmap,save_utmap,save_dudx,save_dudy,save_dudz
  integer :: save_dvdx,save_dvdy,save_dvdz,save_dwdx,save_dwdy,save_dwdz
  integer :: save_dphidx,save_dphidy,save_dphidz,save_abs,save_V

  !module tripping
  integer ::  z_modes, nxt_itr, itrip
  real(mytype) :: x0_tr, xs_tr, ys_tr, ts_tr, zs_param, zs_tr, randomseed, A_trip
  real(mytype), allocatable, dimension(:) :: h_coeff, h_nxt,h_i

  !numbers
  
  real(mytype),parameter :: zpone=0.1_mytype
  real(mytype),parameter :: zptwo=0.2_mytype
  real(mytype),parameter :: zpthree=0.3_mytype
  real(mytype),parameter :: zpfour=0.4_mytype
  real(mytype),parameter :: zpfive=0.5_mytype
  real(mytype),parameter :: zpsix=0.6_mytype
  real(mytype),parameter :: zpseven=0.7_mytype
  real(mytype),parameter :: zpeight=0.8_mytype
  real(mytype),parameter :: zpnine=0.9_mytype

  real(mytype),parameter :: half=0.5_mytype
  real(mytype),parameter :: zero=0._mytype
  real(mytype),parameter :: one=1._mytype
  real(mytype),parameter :: onepfive=1.5_mytype
  real(mytype),parameter :: two=2._mytype
  real(mytype),parameter :: three=3._mytype
  real(mytype),parameter :: four=4._mytype
  real(mytype),parameter :: five=5._mytype
  real(mytype),parameter :: six=6._mytype
  real(mytype),parameter :: seven=7._mytype
  real(mytype),parameter :: eight=8._mytype
  real(mytype),parameter :: nine=9._mytype

  real(mytype),parameter :: ten=10._mytype
  real(mytype),parameter :: eleven=11._mytype
  real(mytype),parameter :: twelve=12._mytype
  real(mytype),parameter :: thirteen=13._mytype
  real(mytype),parameter :: fourteen=14._mytype
  real(mytype),parameter :: fifteen=15._mytype
  real(mytype),parameter :: sixteen=16._mytype
  real(mytype),parameter :: seventeen=17._mytype

  real(mytype),parameter :: twenty=20._mytype
  real(mytype),parameter :: twentyfour=24._mytype
  real(mytype),parameter :: twentyfive=25._mytype
  real(mytype),parameter :: twentyseven=27._mytype
  real(mytype),parameter :: twentyeight=28._mytype
!
  real(mytype),parameter :: thirytwo=32._mytype
  real(mytype),parameter :: thirtysix=36._mytype
!
  real(mytype),parameter :: fortyfour=44._mytype
  real(mytype),parameter :: fortyfive=45._mytype
  real(mytype),parameter :: fortyeight=48._mytype
!
  real(mytype),parameter :: sixty=60._mytype
  real(mytype),parameter :: sixtytwo=62._mytype
  real(mytype),parameter :: sixtythree=63._mytype
!
  real(mytype),parameter :: seventy=70._mytype
  real(mytype),parameter :: seventyfive=75._mytype
!
  real(mytype),parameter :: onehundredtwentysix=126._mytype
  real(mytype),parameter :: onehundredtwentyeight=128._mytype
!
  real(mytype),parameter :: twohundredsix=206._mytype
  real(mytype),parameter :: twohundredeight=208._mytype
  real(mytype),parameter :: twohundredfiftysix=256._mytype
  real(mytype),parameter :: twohundredseventytwo=272._mytype
!
  real(mytype),parameter :: twothousand=2000._mytype
  real(mytype),parameter :: thirtysixthousand=3600._mytype


#ifdef DOUBLE_PREC 
  real(mytype),parameter :: pi=dacos(-one) 
  real(mytype),parameter :: twopi=two*dacos(-one) 
#else
  real(mytype),parameter :: pi=acos(-one)
  real(mytype),parameter :: twopi=two*acos(-one) 
#endif

end module param

module complex_geometry

use decomp_2d,only : mytype
use variables,only : nx,ny,nz,nxm,nym,nzm

  integer     ,allocatable,dimension(:,:)   :: nobjx,nobjy,nobjz
  integer     ,allocatable,dimension(:,:,:) :: nxipif,nxfpif,nyipif,nyfpif,nzipif,nzfpif
  real(mytype),allocatable,dimension(:,:,:) :: xi,xf,yi,yf,zi,zf
  integer :: nxraf,nyraf,nzraf,nraf,nobjmax
end module complex_geometry

module derivX

  use decomp_2d, only : mytype

  real(mytype) :: alcaix6,acix6,bcix6
  real(mytype) :: ailcaix6,aicix6,bicix6,cicix6,dicix6
  real(mytype) :: alfa1x,af1x,bf1x,cf1x,df1x,alfa2x,af2x,alfanx,afnx,bfnx
  real(mytype) :: cfnx,dfnx,alfamx,afmx,alfaix,afix,bfix,alsa1x,as1x,bs1x
  real(mytype) :: cs1x,ds1x,alsa2x,as2x,alsanx,asnx,bsnx,csnx,dsnx,alsamx
  real(mytype) :: asmx,alsa3x,as3x,bs3x,alsatx,astx,bstx

  !O6SVV
  real(mytype) :: alsa4x,as4x,bs4x,cs4x
  real(mytype) :: alsattx,asttx,bsttx,csttx
  real(mytype) :: alsaix,asix,bsix,csix,dsix

end module derivX

module derivY

  use decomp_2d, only : mytype

  real(mytype) :: alcaiy6,aciy6,bciy6
  real(mytype) :: ailcaiy6,aiciy6,biciy6,ciciy6,diciy6
  real(mytype) :: alfa1y,af1y,bf1y,cf1y,df1y,alfa2y,af2y,alfany,afny,bfny
  real(mytype) :: cfny,dfny,alfamy,afmy,alfajy,afjy,bfjy,alsa1y,as1y,bs1y
  real(mytype) :: cs1y,ds1y,alsa2y,as2y,alsany,asny,bsny,csny,dsny,alsamy
  real(mytype) :: asmy,alsa3y,as3y,bs3y,alsaty,asty,bsty

  !O6SVV
  real(mytype) :: alsa4y,as4y,bs4y,cs4y
  real(mytype) :: alsatty,astty,bstty,cstty
  real(mytype) :: alsajy,asjy,bsjy,csjy,dsjy

end module derivY

module derivZ

  use decomp_2d, only : mytype

  real(mytype) :: alcaiz6,aciz6,bciz6
  real(mytype) :: ailcaiz6,aiciz6,biciz6,ciciz6,diciz6
  real(mytype) :: alfa1z,af1z,bf1z,cf1z,df1z,alfa2z,af2z,alfanz,afnz,bfnz
  real(mytype) :: cfnz,dfnz,alfamz,afmz,alfakz,afkz,bfkz,alsa1z,as1z,bs1z
  real(mytype) :: cs1z,ds1z,alsa2z,as2z,alsanz,asnz,bsnz,csnz,dsnz,alsamz
  real(mytype) :: asmz,alsa3z,as3z,bs3z,alsatz,astz,bstz

  !O6SVV
  real(mytype) :: alsa4z,as4z,bs4z,cs4z
  real(mytype) :: alsattz,asttz,bsttz,csttz
  real(mytype) :: alsakz,askz,bskz,cskz,dskz

end module derivZ

module parfiX

  use decomp_2d, only : mytype

  real(mytype) :: fia1x, fib1x, fic1x, fid1x, fie1x, fia2x, fib2x, fic2x, fid2x
  real(mytype) :: fie2x, fia3x, fib3x, fic3x, fid3x, fie3x, fianx, fibnx, ficnx, fidnx
  real(mytype) :: fienx, fiamx, fibmx, ficmx, fidmx, fiemx, fiapx, fibpx, ficpx, fidpx
  real(mytype) :: fiepx, fiaix, fibix, ficix, fidix, fialx, fibex, fih1x, fih2x, fih3x,fih4x
end module parfiX

module parfiY

  use decomp_2d, only : mytype

  real(mytype) :: fia1y, fib1y, fic1y, fid1y, fie1y, fia2y, fib2y, fic2y, fid2y
  real(mytype) :: fie2y, fia3y, fib3y, fic3y, fid3y, fie3y, fiany, fibny, ficny, fidny
  real(mytype) :: fieny, fiamy, fibmy, ficmy, fidmy, fiemy, fiapy, fibpy, ficpy, fidpy
  real(mytype) :: fiepy, fiaiy, fibiy, ficiy, fidiy, fialy, fibey, fih1y, fih2y, fih3y,fih4y
end module parfiY

module parfiZ

  use decomp_2d, only : mytype

  real(mytype) :: fia1z, fib1z, fic1z, fid1z, fie1z, fia2z, fib2z, fic2z, fid2z
  real(mytype) :: fie2z, fia3z, fib3z, fic3z, fid3z, fie3z, fianz, fibnz, ficnz, fidnz
  real(mytype) :: fienz, fiamz, fibmz, ficmz, fidmz, fiemz, fiapz, fibpz, ficpz, fidpz
  real(mytype) :: fiepz, fiaiz, fibiz, ficiz, fidiz, fialz, fibez, fih1z, fih2z, fih3z,fih4z

end module parfiZ
