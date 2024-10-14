#!/bin/bash
if [[ $target == *mingw* ]] ; then

# dopri
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5.o ./dopri5.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5_i32.o ./dopri5.f
gfortran -shared -o ./dopri5.dll ./dopri5.o
gfortran -shared -o ./dopri5_i32.dll ./dopri5_i32.o
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853.o ./dop853.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853_i32.o ./dop853.f
gfortran -shared -o ./dop853.dll ./dop853.o
gfortran -shared -o ./dop853_i32.dll ./dop853_i32.o

# odex
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./odex.o ./odex.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./odex_i32.o ./odex.f
gfortran -shared -o ./odex.dll ./odex.o
gfortran -shared -o ./odex_i32.dll ./odex_i32.o

# lapack
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack.o ./dc_lapack.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack_i32.o ./dc_lapack.f
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack.o ./lapack.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack_i32.o ./lapack.f
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc.o ./lapackc.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc_i32.o ./lapackc.f

# radau
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau5.o ./radau5.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau5_i32.o ./radau5.f
gfortran -shared -o ./radau5.dll ./radau5.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -o ./radau5_i32.dll ./radau5_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau.o ./radau.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau_i32.o ./radau.f
gfortran -shared -o ./radau.dll ./radau.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -o ./radau_i32.dll ./radau_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# seulex
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./seulex.o ./seulex.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./seulex_i32.o ./seulex.f
gfortran -shared -o ./seulex.dll ./seulex.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -o ./seulex_i32.dll ./seulex_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# rodas
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./rodas.o ./rodas.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./rodas_i32.o ./rodas.f
gfortran -shared -o ./rodas.dll ./rodas.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -o ./rodas_i32.dll ./rodas_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# slatec
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec.o ./slatec.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec_i32.o ./slatec.f

# ddeabm
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddeabm.o ./ddeabm.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddeabm_i32.o ./ddeabm.f
gfortran -shared -o ./ddeabm.dll ./ddeabm.o ./slatec.o
gfortran -shared -o ./ddeabm_i32.dll ./ddeabm_i32.o ./slatec_i32.o

# ddebdf
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddebdf.o ./ddebdf.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddebdf_i32.o ./ddebdf.f
gfortran -shared -o ./ddebdf.dll ./ddebdf.o ./slatec.o
gfortran -shared -o ./ddebdf_i32.dll ./ddebdf_i32.o ./slatec_i32.o

# bvpsol
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol.o ./bvpsol.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol_i32.o ./bvpsol.f
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol.o ./linalg_bvpsol.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol_i32.o ./linalg_bvpsol.f
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst.o ./zibconst.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst_i32.o ./zibconst.f
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol.o ./ma28_bvpsol.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol_i32.o ./ma28_bvpsol.f
gfortran -shared -o ./bvpsol.dll ./bvpsol.o ./linalg_bvpsol.o ./zibconst.o ./ma28_bvpsol.o
gfortran -shared -o ./bvpsol_i32.dll ./bvpsol_i32.o ./linalg_bvpsol_i32.o ./zibconst_i32.o ./ma28_bvpsol_i32.o

# colnew
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew.o ./colnew.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew_i32.o ./colnew.f
gfortran -shared -o ./colnew.dll ./colnew.o
gfortran -shared -o ./colnew_i32.dll ./colnew_i32.o

# coldae
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae.o ./coldae.f
gfortran -c -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae_i32.o ./coldae.f
gfortran -shared -o ./coldae.dll ./coldae.o
gfortran -shared -o ./coldae_i32.dll ./coldae_i32.o

# bvpm2
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvp_la-2.o ./bvp_la-2.f
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -std=f2008 -o ./bvp_m-2.o ./bvp_m-2.f90
gfortran -c -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wall -Wextra -Wimplicit-interface -std=f2008ts -o ./bvp_m_proxy.o ./bvp_m_proxy.f90
gfortran -shared -o ./bvp_m_proxy.dll ./bvp_m_proxy.o ./bvp_m-2.o ./bvp_la-2.o
cp *.dll $libdir



elif [[ $target == *apple* ]] ; then

# dopri
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5.o ./dopri5.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5_i32.o ./dopri5.f
gfortran -shared -fPIC -o ./dopri5.dylib ./dopri5.o
gfortran -shared -fPIC -o ./dopri5_i32.dylib ./dopri5_i32.o
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853.o ./dop853.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853_i32.o ./dop853.f
gfortran -shared -fPIC -o ./dop853.dylib ./dop853.o
gfortran -shared -fPIC -o ./dop853_i32.dylib ./dop853_i32.o

# odex
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./odex.o ./odex.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./odex_i32.o ./odex.f
gfortran -shared -fPIC -o ./odex.dylib ./odex.o
gfortran -shared -fPIC -o ./odex_i32.dylib ./odex_i32.o

# lapack
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack.o ./dc_lapack.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack_i32.o ./dc_lapack.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack.o ./lapack.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack_i32.o ./lapack.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc.o ./lapackc.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc_i32.o ./lapackc.f

# radau
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau5.o ./radau5.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau5_i32.o ./radau5.f
gfortran -shared -fPIC -o ./radau5.dylib ./radau5.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./radau5_i32.dylib ./radau5_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau.o ./radau.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau_i32.o ./radau.f
gfortran -shared -fPIC -o ./radau.dylib ./radau.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./radau_i32.dylib ./radau_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# seulex
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./seulex.o ./seulex.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./seulex_i32.o ./seulex.f
gfortran -shared -fPIC -o ./seulex.dylib ./seulex.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./seulex_i32.dylib ./seulex_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# rodas
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./rodas.o ./rodas.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./rodas_i32.o ./rodas.f
gfortran -shared -fPIC -o ./rodas.dylib ./rodas.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./rodas_i32.dylib ./rodas_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# slatec
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec.o ./slatec.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec_i32.o ./slatec.f

# ddeabm
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddeabm.o ./ddeabm.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddeabm_i32.o ./ddeabm.f
gfortran -shared -fPIC -o ./ddeabm.dylib ./ddeabm.o ./slatec.o
gfortran -shared -fPIC -o ./ddeabm_i32.dylib ./ddeabm_i32.o ./slatec_i32.o

# ddebdf
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddebdf.o ./ddebdf.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddebdf_i32.o ./ddebdf.f
gfortran -shared -fPIC -o ./ddebdf.dylib ./ddebdf.o ./slatec.o
gfortran -shared -fPIC -o ./ddebdf_i32.dylib ./ddebdf_i32.o ./slatec_i32.o

# bvpsol
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol.o ./bvpsol.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol_i32.o ./bvpsol.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol.o ./linalg_bvpsol.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol_i32.o ./linalg_bvpsol.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst.o ./zibconst.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst_i32.o ./zibconst.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol.o ./ma28_bvpsol.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol_i32.o ./ma28_bvpsol.f
gfortran -shared -fPIC -o ./bvpsol.dylib ./bvpsol.o ./linalg_bvpsol.o ./zibconst.o ./ma28_bvpsol.o
gfortran -shared -fPIC -o ./bvpsol_i32.dylib ./bvpsol_i32.o ./linalg_bvpsol_i32.o ./zibconst_i32.o ./ma28_bvpsol_i32.o

# colnew
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew.o ./colnew.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew_i32.o ./colnew.f
gfortran -shared -fPIC -o ./colnew.dylib ./colnew.o
gfortran -shared -fPIC -o ./colnew_i32.dylib ./colnew_i32.o

# coldae
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae.o ./coldae.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae_i32.o ./coldae.f
gfortran -shared -fPIC -o ./coldae.dylib ./coldae.o
gfortran -shared -fPIC -o ./coldae_i32.dylib ./coldae_i32.o

# bvpm2
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvp_la-2.o ./bvp_la-2.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -std=f2008 -o ./bvp_m-2.o ./bvp_m-2.f90
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wall -Wextra -Wimplicit-interface -std=f2008ts -o ./bvp_m_proxy.o ./bvp_m_proxy.f90
gfortran -shared -fPIC -o ./bvp_m_proxy.dylib ./bvp_m_proxy.o ./bvp_m-2.o ./bvp_la-2.o
cp *.dylib $libdir



else

# dopri
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5.o ./dopri5.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dopri5_i32.o ./dopri5.f
gfortran -shared -fPIC -o ./dopri5.so ./dopri5.o
gfortran -shared -fPIC -o ./dopri5_i32.so ./dopri5_i32.o
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853.o ./dop853.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dop853_i32.o ./dop853.f
gfortran -shared -fPIC -o ./dop853.so ./dop853.o
gfortran -shared -fPIC -o ./dop853_i32.so ./dop853_i32.o

# odex
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./odex.o ./odex.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./odex_i32.o ./odex.f
gfortran -shared -fPIC -o ./odex.so ./odex.o
gfortran -shared -fPIC -o ./odex_i32.so ./odex_i32.o

# lapack
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack.o ./dc_lapack.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./dc_lapack_i32.o ./dc_lapack.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack.o ./lapack.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapack_i32.o ./lapack.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc.o ./lapackc.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./lapackc_i32.o ./lapackc.f

# radau
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau5.o ./radau5.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau5_i32.o ./radau5.f
gfortran -shared -fPIC -o ./radau5.so ./radau5.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./radau5_i32.so ./radau5_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau.o ./radau.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./radau_i32.o ./radau.f
gfortran -shared -fPIC -o ./radau.so ./radau.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./radau_i32.so ./radau_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# seulex
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./seulex.o ./seulex.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./seulex_i32.o ./seulex.f
gfortran -shared -fPIC -o ./seulex.so ./seulex.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./seulex_i32.so ./seulex_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# rodas
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./rodas.o ./rodas.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./rodas_i32.o ./rodas.f
gfortran -shared -fPIC -o ./rodas.so ./rodas.o ./dc_lapack.o ./lapack.o ./lapackc.o
gfortran -shared -fPIC -o ./rodas_i32.so ./rodas_i32.o ./dc_lapack_i32.o ./lapack_i32.o ./lapackc_i32.o

# slatec
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec.o ./slatec.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./slatec_i32.o ./slatec.f

# ddeabm
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddeabm.o ./ddeabm.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddeabm_i32.o ./ddeabm.f
gfortran -shared -fPIC -o ./ddeabm.so ./ddeabm.o ./slatec.o
gfortran -shared -fPIC -o ./ddeabm_i32.so ./ddeabm_i32.o ./slatec_i32.o

# ddebdf
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddebdf.o ./ddebdf.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ddebdf_i32.o ./ddebdf.f
gfortran -shared -fPIC -o ./ddebdf.so ./ddebdf.o ./slatec.o
gfortran -shared -fPIC -o ./ddebdf_i32.so ./ddebdf_i32.o ./slatec_i32.o

# bvpsol
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol.o ./bvpsol.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvpsol_i32.o ./bvpsol.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol.o ./linalg_bvpsol.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./linalg_bvpsol_i32.o ./linalg_bvpsol.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst.o ./zibconst.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./zibconst_i32.o ./zibconst.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol.o ./ma28_bvpsol.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./ma28_bvpsol_i32.o ./ma28_bvpsol.f
gfortran -shared -fPIC -o ./bvpsol.so ./bvpsol.o ./linalg_bvpsol.o ./zibconst.o ./ma28_bvpsol.o
gfortran -shared -fPIC -o ./bvpsol_i32.so ./bvpsol_i32.o ./linalg_bvpsol_i32.o ./zibconst_i32.o ./ma28_bvpsol_i32.o

# colnew
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew.o ./colnew.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./colnew_i32.o ./colnew.f
gfortran -shared -fPIC -o ./colnew.so ./colnew.o
gfortran -shared -fPIC -o ./colnew_i32.so ./colnew_i32.o

# coldae
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae.o ./coldae.f
gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./coldae_i32.o ./coldae.f
gfortran -shared -fPIC -o ./coldae.so ./coldae.o
gfortran -shared -fPIC -o ./coldae_i32.so ./coldae_i32.o

# bvpm2
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -o ./bvp_la-2.o ./bvp_la-2.f
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -std=f2008 -o ./bvp_m-2.o ./bvp_m-2.f90
gfortran -c -fPIC -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -Wall -Wextra -Wimplicit-interface -std=f2008ts -o ./bvp_m_proxy.o ./bvp_m_proxy.f90
gfortran -shared -fPIC -o ./bvp_m_proxy.so ./bvp_m_proxy.o ./bvp_m-2.o ./bvp_la-2.o
cp *.so $libdir



fi

