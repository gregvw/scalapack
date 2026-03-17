add_subdirectory(LAPACK)

set (ATOOLS
         iceil.f      ilacpy.f     ilcm.f       indxg2p.f    indxg2l.f
         indxl2g.f    infog1l.f    infog2l.f    npreroc.f    numroc.f
         chk1mat.f    pchkxmat.f   sltimer.f    desc_convert.f
         descinit.f   descset.f
         SL_init.f)

set (ATOOLS_I8
         numroc_i8.f    descinit_i8.f  descset_i8.f   infog2l_i8.f
         infog1l_i8.f   indxl2g_i8.f   indxg2p_i8.f   indxg2l_i8.f
         chk1mat_i8.f   pchkxmat_i8.f  desc_convert_i8.f
         pdelset_i8.f   pselset_i8.f   pcelset_i8.f   pzelset_i8.f)

set (ITOOLS
         picol2row.f  pirow2col.f  pilaprnt.f   pitreecomb.f  pifillpad.f 
         pichekpad.f  pielset.f    pielset2.f   pielget.f)

set (STOOLS
         dsnrm2.f     dsasum.f     slatcpy.f    ssdot.f       smatadd.f  
         psmatadd.f   pscol2row.f  psrow2col.f  pslaprnt.f    pstreecomb.f
         psfillpad.f  pschekpad.f  pselset.f    pselset2.f    pselget.f
         pslaread.f   pslawrite.f)

set (DTOOLS
         dddot.f      dlatcpy.f    dmatadd.f    pdmatadd.f    pdcol2row.f
         pdrow2col.f  pdlaprnt.f   pdtreecomb.f pdfillpad.f   pdchekpad.f
         pdelset.f    pdelset2.f   pdelget.f
         pdlaread.f   pdlawrite.f)

set (CTOOLS
         dscnrm2.f    dscasum.f    ccdotu.f     ccdotc.f      clatcpy.f
         cmatadd.f    pcmatadd.f   pccol2row.f  pcrow2col.f   pclaprnt.f
         pctreecomb.f pcfillpad.f  pcchekpad.f  pcelset.f     pcelset2.f
         pcelget.f 
         pclaread.f   pclawrite.f)

set (ZTOOLS
         zzdotu.f     zzdotc.f     zlatcpy.f    zmatadd.f     pzmatadd.f
         pzcol2row.f  pzrow2col.f  pzlaprnt.f   pztreecomb.f  pzfillpad.f
         pzchekpad.f  pzelset.f    pzelset2.f   pzelget.f
         pzlaread.f   pzlawrite.f)

set(tools
     ${ATOOLS} ${ATOOLS_I8} ${ITOOLS} ${STOOLS} ${DTOOLS} ${CTOOLS} ${ZTOOLS})
     
set(tools-C
    reshape.c  SL_gridreshape.c  globchk_i8.c )

# Standalone I8 tools test (no MPI required)
add_executable(xi8tools test_i8_tools.f)
target_link_libraries(xi8tools scalapack)
add_test(NAME xi8tools COMMAND xi8tools)
     
