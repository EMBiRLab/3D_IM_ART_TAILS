function ptau_plen = two_link_2JDOF_varied_len_ptau_plen(in1,in2,in3,in4)
%two_link_2JDOF_varied_len_ptau_plen
%    PTAU_PLEN = two_link_2JDOF_varied_len_ptau_plen(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    24-Jan-2024 11:18:55

l1 = in4(1,:);
l2 = in4(2,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
qd1 = in2(1,:);
qd2 = in2(2,:);
qd3 = in2(3,:);
qd4 = in2(4,:);
qd5 = in2(5,:);
qd6 = in2(6,:);
qd7 = in2(7,:);
qdd1 = in3(1,:);
qdd2 = in3(2,:);
qdd3 = in3(3,:);
qdd4 = in3(4,:);
qdd5 = in3(5,:);
qdd6 = in3(6,:);
qdd7 = in3(7,:);
t2 = cos(q2);
t3 = cos(q3);
t4 = cos(q4);
t5 = cos(q5);
t6 = cos(q6);
t7 = cos(q7);
t8 = sin(q2);
t9 = sin(q3);
t10 = sin(q4);
t11 = sin(q5);
t12 = sin(q6);
t13 = sin(q7);
t14 = l1.^2;
t15 = l2.^2;
t27 = qd3./2.0;
t28 = qdd3./2.0;
t16 = qd2.*t3;
t17 = qdd1.*t2;
t18 = qdd2.*t3;
t19 = qd1.*t8;
t20 = qd2.*t9;
t21 = qdd1.*t8;
t22 = qdd2.*t9;
t23 = qd1.*qd2.*t2;
t25 = t14.*2.0;
t26 = t15.*2.0;
t29 = qd1.*t2.*t3;
t31 = qd1.*t2.*t9;
t24 = qd2.*t19;
t30 = qd3+t19;
t32 = -t18;
t34 = t19./2.0;
t35 = t21./2.0;
t36 = t23./2.0;
t37 = -t31;
t40 = t20+t29;
t41 = qdd3+t21+t23;
t43 = t25+1.0./6.0e+2;
t44 = t26+1.0./6.0e+2;
t33 = -t24;
t38 = t4.*t30;
t39 = t10.*t30;
t45 = qd4+t40;
t46 = qd3.*t40;
t47 = t16+t37;
t48 = t27+t34;
t49 = t4.*t40;
t50 = t4.*t41;
t51 = t10.*t41;
t54 = qd4.*t10.*t40;
t62 = t10.*t11.*t40;
t74 = (t5.*t10.*t40)./2.0;
t78 = t28+t35+t36;
t42 = t17+t33;
t52 = l1.*t49;
t53 = t5.*t45;
t55 = t11.*t45;
t56 = qd3.*t47;
t59 = t5.*t48;
t60 = t11.*t48;
t61 = t6.*t49;
t63 = t4.*t47;
t64 = t10.*t47;
t66 = -t54;
t68 = (qd4.*t49)./2.0;
t69 = t54./2.0;
t70 = -t62;
t75 = (t12.*t49)./2.0;
t76 = t62./2.0;
t79 = t5.*t78;
t80 = t11.*t78;
t57 = t3.*t42;
t58 = t9.*t42;
t65 = t59.*2.0;
t67 = -t55;
t71 = -t64;
t72 = -t69;
t73 = t61./2.0;
t77 = -t76;
t81 = t79.*2.0;
t82 = t80.*2.0;
t83 = -t79;
t86 = t39+t63;
t105 = t60+t74;
t84 = -t81;
t85 = -t82;
t87 = qd4.*t86;
t88 = t38+t71;
t89 = t5.*t86;
t90 = t11.*t86;
t103 = t32+t46+t58;
t104 = t22+t56+t57;
t107 = t59+t77;
t108 = qd5.*t105;
t114 = t6.*t105;
t115 = t12.*t105;
t91 = qd5+t88;
t92 = qd4.*t88;
t93 = -t87;
t106 = qdd4+t104;
t109 = t4.*t103;
t110 = t4.*t104;
t111 = t10.*t103;
t112 = qd5.*t107;
t113 = t108.*2.0;
t120 = t115.*2.0;
t121 = l1.*t107.*2.0;
t124 = (t10.*t104)./2.0;
t125 = -t115;
t129 = t53+t90;
t133 = t67+t89;
t141 = -qd5.*(t55-t89);
t142 = -t12.*(t55-t89);
t147 = -t6.*(t55-t89);
t151 = t12.*(t55-t89);
t157 = -t49.*(t55-t89);
t159 = -t52.*(t55-t89);
t202 = t107.*(t55-t89).*-2.0;
t94 = l2.*t91;
t95 = t91.^2;
t96 = l1.*t91.*2.0;
t97 = t6.*t91;
t98 = t7.*t91;
t99 = t12.*t91;
t100 = t13.*t91;
t116 = t5.*t106;
t117 = t11.*t106;
t118 = t112.*2.0;
t119 = -t109;
t123 = t110./2.0;
t126 = -t120;
t127 = t43.*t91;
t130 = qd6+t129;
t131 = t129.^2;
t132 = qd5.*t129;
t134 = l1.*t129.*2.0;
t135 = l2.*t6.*t129;
t136 = qd6.*t6.*t129;
t138 = l2.*t12.*t129;
t139 = qd6.*t12.*t129;
t144 = t7.*t12.*t129;
t145 = t12.*t13.*t129;
t150 = t6.*t15.*t129;
t154 = t43.*t129;
t155 = t49.*t129;
t156 = t68+t124;
t171 = t91.*t105.*2.0;
t175 = t91.*t107.*2.0;
t177 = l1.*t91.*t105.*-2.0;
t183 = l1.*t91.*(t55-t89).*-2.0;
t184 = t91.*(t55-t89).*(-1.0./3.0e+2);
t185 = t105.*t129.*2.0;
t186 = (t91.*(t55-t89))./3.0e+2;
t190 = qdd5+t50+t93+t111;
t210 = l1.*t202;
t230 = l1.*t129.*(t55-t89).*-2.0;
t232 = t129.*(t55-t89).*(-1.0./3.0e+2);
t235 = (t129.*(t55-t89))./3.0e+2;
t101 = l1.*t95.*2.0;
t102 = -t97;
t122 = -t118;
t128 = t117./3.0e+2;
t137 = t7.*t130;
t140 = t13.*t130;
t143 = t135.*4.0;
t146 = l1.*t131.*2.0;
t148 = -t139;
t152 = t132./3.0e+2;
t153 = -t145;
t158 = t72+t123;
t160 = t49+t134;
t161 = t5.*t156;
t162 = t11.*t156;
t163 = l1.*(t69-t123).*-2.0;
t164 = -t6.*(t69-t123);
t166 = -t12.*(t69-t123);
t169 = t6.*(t69-t123).*-2.0;
t170 = t6.*(t69-t123);
t172 = t94+t107;
t173 = t65+t70+t96;
t174 = t96.*t105;
t176 = -t171;
t182 = t52+t154;
t187 = t121+t127;
t188 = t100+t144;
t189 = t99+t147;
t191 = t105.*t134;
t192 = t51+t92+t119;
t201 = l2.*t190;
t203 = qd7+t97+t151;
t204 = l1.*t190.*2.0;
t205 = t6.*t190;
t206 = t7.*t190;
t208 = t12.*t190;
t209 = t13.*t190;
t211 = -qd6.*(t97+t151);
t220 = -t127.*(t55-t89);
t231 = t43.*t190;
t238 = -t154.*(t55-t89);
t239 = t75+t114+t138;
t241 = t73+t125+t135;
t149 = -t140;
t165 = t161.*2.0;
t167 = t162.*2.0;
t178 = t7.*t172;
t179 = t13.*t172;
t193 = (t49.*t173)./2.0;
t194 = qd7.*t188;
t195 = qd6.*t189;
t196 = t98+t153;
t197 = t102+t142;
t199 = t7.*t189;
t200 = t13.*t189;
t213 = t5.*t192;
t214 = -t204;
t215 = -t201;
t216 = t11.*t192;
t218 = t203.^2;
t221 = l2.*t203.*2.0;
t222 = -t206;
t223 = -t208;
t225 = t15.*t203;
t233 = t44.*t203;
t234 = -t231;
t236 = t107.*t160;
t237 = t91.*t182;
t240 = qd6.*t239;
t242 = t7.*t239;
t243 = t13.*t239;
t246 = qd6.*t241;
t249 = l2.*t241.*2.0;
t250 = t80+t112+t161;
t251 = t83+t108+t162;
t252 = t129.*t187;
t270 = l2.*t188.*t203.*-2.0;
t301 = -t13.*(t79-t108-t162+t201);
t302 = -t7.*(t79-t108-t162+t201);
t303 = t13.*(t79-t108-t162+t201).*-2.0;
t304 = t7.*(t79-t108-t162+t201).*-2.0;
t168 = -t165;
t180 = t178.*2.0;
t181 = -t178;
t198 = -t193;
t207 = qd7.*t196;
t212 = -t195;
t217 = l2.*t196.*2.0;
t219 = t15.*t196;
t224 = -t213;
t226 = l2.*t218.*2.0;
t227 = t15.*t218;
t228 = t213./3.0e+2;
t244 = t240.*2.0;
t245 = t243.*2.0;
t253 = t6.*t250;
t254 = t12.*t250;
t255 = l1.*t251.*2.0;
t257 = -t252;
t259 = t137+t200;
t260 = t149+t199;
t264 = -qd7.*(t140-t199);
t269 = t188.*t221;
t271 = t188.*t225;
t274 = t116+t141+t216;
t286 = t179+t242;
t293 = -qd7.*(t178-t243);
t294 = t135.*(t140-t199).*-2.0;
t295 = l2.*(t178-t243).*-2.0;
t296 = qd7.*(t178-t243).*-2.0;
t297 = -t150.*(t140-t199);
t298 = l2.*(t178-t243).*2.0;
t300 = t215+t251;
t307 = l2.*t196.*(t140-t199).*-2.0;
t309 = l2.*t203.*(t140-t199).*-2.0;
t310 = -t225.*(t140-t199);
t312 = t203.*(t140-t199).*(-1.0./3.0e+2);
t313 = (t203.*(t140-t199))./3.0e+2;
t319 = -t233.*(t140-t199);
t323 = t241.*(t140-t199).*-2.0;
t330 = t203.*(t178-t243).*-2.0;
t332 = t203.*(t178-t243).*2.0;
t333 = t221.*(t178-t243);
t343 = (t140-t199).*(t178-t243).*2.0;
t375 = t84+t113+t157+t167+t176+t214+t230;
t376 = -t11.*(t81-t113-t167+t171+t204+t49.*(t55-t89)+t134.*(t55-t89));
t377 = t11.*(t81-t113-t167+t171+t204+t49.*(t55-t89)+t134.*(t55-t89));
t378 = t5.*(t81-t113-t167+t171+t204+t49.*(t55-t89)+t134.*(t55-t89)).*(-1.0./2.0);
t229 = -t228;
t247 = -t244;
t248 = -t245;
t256 = t254.*2.0;
t261 = t259.^2;
t262 = l2.*t259.*2.0;
t263 = t15.*t259;
t267 = (qd7.*t259)./3.0e+2;
t268 = t44.*t259;
t272 = t203.*t217;
t273 = -t271;
t275 = t117+t132+t224;
t276 = qdd6+t274;
t277 = l1.*t274.*2.0;
t278 = t6.*t274;
t280 = t12.*t274;
t288 = t181+t243;
t289 = qd7.*t286;
t290 = t135.*t259.*2.0;
t299 = t43.*t274;
t308 = -t219.*(t140-t199);
t321 = t241.*t259.*2.0;
t322 = t249.*t259;
t324 = l2.*t323;
t326 = t203.*t286.*2.0;
t327 = t221.*t286;
t329 = l2.*t203.*t286.*-2.0;
t331 = t203.*t295;
t334 = l2.*t259.*(t140-t199).*-2.0;
t336 = t259.*(t140-t199).*(-1.0./3.0e+2);
t337 = t225+t298;
t338 = (t259.*(t140-t199))./3.0e+2;
t341 = t259.*t286.*2.0;
t344 = -t343;
t345 = t298.*(t140-t199);
t346 = t295.*(t140-t199);
t350 = t85+t101+t122+t146+t155+t168+t175;
t365 = t219+t233+t298;
t395 = t159+t177+t234+t235+t238+t255;
t396 = -t4.*(t174+t231+t232-t255+t52.*(t55-t89)+t154.*(t55-t89));
t397 = -t10.*(t174+t231+t232-t255+t52.*(t55-t89)+t154.*(t55-t89));
t258 = -t256;
t265 = l2.*t261.*2.0;
t266 = t15.*t261;
t279 = t6.*t275;
t281 = t12.*t275;
t282 = l2.*t278;
t283 = t7.*t276;
t284 = l2.*t280;
t291 = (t13.*t276)./3.0e+2;
t292 = t289.*2.0;
t305 = t188.*t262;
t306 = t188.*t263;
t311 = t136+t280;
t314 = t148+t278;
t317 = l2.*(t139-t278).*-2.0;
t318 = -t15.*(t139-t278);
t320 = t61+t126+t143+t262;
t325 = t249+t263;
t328 = -t326;
t335 = -t263.*(t140-t199);
t339 = -t268.*(t140-t199);
t340 = t150+t249+t268;
t342 = t262.*t286;
t347 = t6.*t129.*t337;
t351 = t5.*t350;
t352 = (t11.*t350)./2.0;
t361 = t180+t217+t221+t248;
t381 = t66+t110+t183+t185+t202+t277;
t383 = t4.*(t54-t110-t185-t277+t96.*(t55-t89)+t107.*(t55-t89).*2.0).*(-1.0./2.0);
t393 = t259.*t365;
t398 = t163+t186+t191+t210+t220+t299;
t399 = -t5.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0);
t400 = -t11.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0);
t401 = t128+t152+t198+t229+t236+t237+t257;
t285 = t282.*2.0;
t287 = -t282;
t315 = t13.*t311;
t316 = t7.*t311;
t348 = -t347;
t349 = t196.*t325;
t353 = qdd7+t205+t212+t281;
t354 = t211+t223+t279;
t357 = -t13.*(t208-t279+qd6.*(t97+t151));
t359 = t13.*(t208-t279+qd6.*(t97+t151));
t362 = t7.*(t208-t279+qd6.*(t97+t151)).*(-1.0./3.0e+2);
t366 = t203.*t340;
t373 = -t320.*(t178-t243);
t374 = t320.*(t178-t243);
t379 = t241.*t361;
t382 = t166+t246+t253+t284;
t394 = -t393;
t402 = t11.*t401;
t403 = t5.*t401;
t409 = t351+t377;
t355 = l2.*t353.*2.0;
t358 = t15.*t353;
t363 = t44.*t353;
t367 = t194+t222+t315;
t368 = t207+t209+t316;
t380 = -t379;
t384 = t170+t240+t254+t287;
t385 = t7.*t382;
t386 = t13.*t382;
t404 = -t403;
t405 = t264+t283+t359;
t410 = (t10.*t409)./2.0;
t417 = -t4.*(t403+t11.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0));
t418 = -t10.*(t403+t11.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0));
t356 = -t355;
t360 = -t358;
t364 = -t363;
t369 = l2.*t367.*2.0;
t370 = l2.*t368.*2.0;
t372 = t15.*t367;
t387 = t385.*2.0;
t388 = t386.*2.0;
t389 = -t385;
t391 = l2.*t384.*2.0;
t406 = l2.*t405.*2.0;
t407 = t15.*t405;
t408 = t44.*t405;
t411 = -t410;
t412 = t289+t302+t386;
t415 = l2.*(t385+qd7.*(t178-t243)+t13.*(t79-t108-t162+t201)).*-2.0;
t416 = t400+t404;
t422 = t397+t417;
t449 = t267+t291+t348+t349+t362+t366+t374+t380+t394;
t371 = -t370;
t390 = -t387;
t392 = -t391;
t413 = l2.*t412.*2.0;
t414 = t293+t301+t389;
t419 = t227+t266+t322+t333+t415;
t424 = -t6.*(t345+t391-t407+t225.*(t140-t199)-l2.*t259.*t286.*2.0);
t432 = t383+t399+t402+t411;
t433 = t169+t247+t258+t285+t305+t307+t309+t317+t341+t344+t406;
t434 = -t6.*(t170.*2.0+t244+t256-t285-t341+t343-t406+l2.*(t139-t278).*2.0+t217.*(t140-t199)+t221.*(t140-t199)-l2.*t188.*t259.*2.0);
t435 = -t12.*(t170.*2.0+t244+t256-t285-t341+t343-t406+l2.*(t139-t278).*2.0+t217.*(t140-t199)+t221.*(t140-t199)-l2.*t188.*t259.*2.0);
t437 = t6.*(t170.*2.0+t244+t256-t285-t341+t343-t406+l2.*(t139-t278).*2.0+t217.*(t140-t199)+t221.*(t140-t199)-l2.*t188.*t259.*2.0);
t438 = t270+t292+t294+t304+t323+t328+t334+t356+t369+t388;
t439 = -t7.*(t269-t292+t326+t355-t369-t388+t135.*(t140-t199).*2.0+t241.*(t140-t199).*2.0+t262.*(t140-t199)+t7.*(t79-t108-t162+t201).*2.0);
t440 = -t13.*(t269-t292+t326+t355-t369-t388+t135.*(t140-t199).*2.0+t241.*(t140-t199).*2.0+t262.*(t140-t199)+t7.*(t79-t108-t162+t201).*2.0);
t441 = t13.*(t269-t292+t326+t355-t369-t388+t135.*(t140-t199).*2.0+t241.*(t140-t199).*2.0+t262.*(t140-t199)+t7.*(t79-t108-t162+t201).*2.0);
t443 = -t7.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0);
t444 = -t13.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0);
t445 = t13.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0);
t450 = t7.*t449;
t451 = t13.*t449;
t420 = t7.*t419;
t421 = t13.*t419;
t423 = t310+t342+t346+t392+t407;
t425 = t324+t329+t335+t360+t413;
t426 = -t7.*(t327+t358-t413+t249.*(t140-t199)+t263.*(t140-t199));
t427 = -t13.*(t327+t358-t413+t249.*(t140-t199)+t263.*(t140-t199));
t428 = t13.*(t327+t358-t413+t249.*(t140-t199)+t263.*(t140-t199));
t429 = t226+t265+t272+t290+t296+t303+t321+t332+t371+t390;
t436 = l2.*t434;
t442 = t306+t308+t313+t318+t319+t342+t346+t392+t408;
t446 = t273+t297+t324+t329+t338+t339+t364+t372+t413;
t447 = -t6.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199));
t448 = -t12.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199));
t464 = t445+t450;
t430 = t7.*t429;
t431 = t13.*t429;
t452 = t420+t428;
t465 = t6.*t464;
t466 = t12.*t464;
t453 = t12.*t452;
t455 = t431+t439;
t456 = t430+t441;
t467 = -t465;
t475 = -t5.*(t465+t12.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199)));
t476 = -t11.*(t465+t12.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199)));
t454 = -t453;
t457 = l2.*t455;
t458 = t11.*t455;
t459 = t6.*t456;
t460 = t12.*t456;
t461 = (t5.*t455)./2.0;
t474 = t448+t467;
t462 = l2.*t460;
t468 = t435+t459;
t469 = t437+t460;
t480 = t421+t426+t447+t457+t466;
t481 = t4.*(t421+t447+t457+t466-t7.*(t327+t358-t413+t249.*(t140-t199)+t263.*(t140-t199)));
t482 = t10.*(t421+t447+t457+t466-t7.*(t327+t358-t413+t249.*(t140-t199)+t263.*(t140-t199)));
t463 = -t462;
t470 = t5.*t468;
t472 = (t11.*t468)./2.0;
t473 = (t4.*t469)./2.0;
t484 = -t5.*(-t451+t453+t462+t6.*(t345+t391-t407+t225.*(t140-t199)-l2.*t259.*t286.*2.0)+l2.*t437+t7.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0));
ptau_plen = ft_1({l1,l2,t10,t105,t107,t11,t110,t113,t12,t121,t123,t127,t129,t134,t139,t140,t15,t150,t154,t167,t171,t174,t184,t185,t199,t2,t204,t219,t225,t231,t232,t233,t249,t255,t259,t263,t268,t271,t277,t278,t286,t299,t3,t306,t312,t327,t336,t345,t352,t358,t363,t372,t378,t391,t396,t399,t4,t402,t403,t407,t408,t410,t413,t418,t421,t424,t436,t437,t443,t451,t453,t454,t457,t458,t461,t462,t463,t465,t466,t470,t472,t473,t475,t476,t481,t482,t484,t49,t5,t52,t54,t55,t6,t69,t7,t8,t81,t89,t9,t96});
end
function ptau_plen = ft_1(ct)
l1 = ct{1};
l2 = ct{2};
t10 = ct{3};
t105 = ct{4};
t107 = ct{5};
t11 = ct{6};
t110 = ct{7};
t113 = ct{8};
t12 = ct{9};
t121 = ct{10};
t123 = ct{11};
t127 = ct{12};
t129 = ct{13};
t134 = ct{14};
t139 = ct{15};
t140 = ct{16};
t15 = ct{17};
t150 = ct{18};
t154 = ct{19};
t167 = ct{20};
t171 = ct{21};
t174 = ct{22};
t184 = ct{23};
t185 = ct{24};
t199 = ct{25};
t2 = ct{26};
t204 = ct{27};
t219 = ct{28};
t225 = ct{29};
t231 = ct{30};
t232 = ct{31};
t233 = ct{32};
t249 = ct{33};
t255 = ct{34};
t259 = ct{35};
t263 = ct{36};
t268 = ct{37};
t271 = ct{38};
t277 = ct{39};
t278 = ct{40};
t286 = ct{41};
t299 = ct{42};
t3 = ct{43};
t306 = ct{44};
t312 = ct{45};
t327 = ct{46};
t336 = ct{47};
t345 = ct{48};
t352 = ct{49};
t358 = ct{50};
t363 = ct{51};
t372 = ct{52};
t378 = ct{53};
t391 = ct{54};
t396 = ct{55};
t399 = ct{56};
t4 = ct{57};
t402 = ct{58};
t403 = ct{59};
t407 = ct{60};
t408 = ct{61};
t410 = ct{62};
t413 = ct{63};
t418 = ct{64};
t421 = ct{65};
t424 = ct{66};
t436 = ct{67};
t437 = ct{68};
t443 = ct{69};
t451 = ct{70};
t453 = ct{71};
t454 = ct{72};
t457 = ct{73};
t458 = ct{74};
t461 = ct{75};
t462 = ct{76};
t463 = ct{77};
t465 = ct{78};
t466 = ct{79};
t470 = ct{80};
t472 = ct{81};
t473 = ct{82};
t475 = ct{83};
t476 = ct{84};
t481 = ct{85};
t482 = ct{86};
t484 = ct{87};
t49 = ct{88};
t5 = ct{89};
t52 = ct{90};
t54 = ct{91};
t55 = ct{92};
t6 = ct{93};
t69 = ct{94};
t7 = ct{95};
t8 = ct{96};
t81 = ct{97};
t89 = ct{98};
t9 = ct{99};
t96 = ct{100};
t485 = -t11.*(-t451+t453+t462+t6.*(t345+t391-t407+t225.*(t140-t199)-l2.*t259.*t286.*2.0)+l2.*t437+t7.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0));
t486 = t5.*(-t451+t453+t462+t6.*(t345+t391-t407+t225.*(t140-t199)-l2.*t259.*t286.*2.0)+l2.*t437+t7.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0));
t488 = -t4.*(t11.*(-t451+t453+t462+t6.*(t345+t391-t407+t225.*(t140-t199)-l2.*t259.*t286.*2.0)+l2.*t437+t7.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0))+t5.*(t465+t12.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199))));
t489 = -t10.*(t11.*(-t451+t453+t462+t6.*(t345+t391-t407+t225.*(t140-t199)-l2.*t259.*t286.*2.0)+l2.*t437+t7.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0))+t5.*(t465+t12.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199))));
t471 = -t470;
t483 = t424+t436+t443+t451+t454+t463;
t487 = t475+t485;
t490 = t482+t488;
t477 = t458+t471;
t478 = (t10.*t477)./2.0;
t479 = -t478;
t491 = t473+t476+t479+t486;
et1 = -t2.*(t3.*(-t402+t410+(t4.*(t54-t110-t185-t277+t96.*(t55-t89)+t107.*(t55-t89).*2.0))./2.0+t5.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0))-t9.*(t4.*(t403+t11.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0))+t10.*(t174+t231+t232-t255+t52.*(t55-t89)+t154.*(t55-t89))));
et2 = -t8.*(t352+t378+t396+t10.*(t403+t11.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0)));
mt1 = [et1+et2,-t9.*(-t402+t410+(t4.*(t54-t110-t185-t277+t96.*(t55-t89)+t107.*(t55-t89).*2.0))./2.0+t5.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0))-t3.*(t4.*(t403+t11.*(t184-t299+l1.*(t69-t123).*2.0+t121.*(t55-t89)+t127.*(t55-t89)-l1.*t105.*t129.*2.0))+t10.*(t174+t231+t232-t255+t52.*(t55-t89)+t154.*(t55-t89)))];
mt2 = [-t352+t418+(t5.*(t81-t113-t167+t171+t204+t49.*(t55-t89)+t134.*(t55-t89)))./2.0+t4.*(t174+t231+t232-t255+t52.*(t55-t89)+t154.*(t55-t89)),t399+t402,t174+t231+t232-t255+t52.*(t55-t89)+t154.*(t55-t89),0.0,0.0];
mt3 = [-t8.*(t461+t472+t481+t10.*(t11.*(-t451+t453+t462+t6.*(t345+t391-t407+t225.*(t140-t199)-l2.*t259.*t286.*2.0)+l2.*t437+t7.*(-t306+t312+t345+t391-t408+t15.*(t139-t278)+t219.*(t140-t199)+t233.*(t140-t199)-l2.*t259.*t286.*2.0))+t5.*(t465+t12.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199)))))-t2.*(t3.*t491+t9.*t490),t3.*t490-t9.*t491,-t461-t472-t481+t489];
mt4 = [t484+t11.*(t465+t12.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199))),-t421-t457-t466+t6.*(t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199))+t7.*(t327+t358-t413+t249.*(t140-t199)+t263.*(t140-t199)),t443+t451,t271+t327+t336+t363-t372-t413+t150.*(t140-t199)+t249.*(t140-t199)+t268.*(t140-t199)];
ptau_plen = reshape([mt1,mt2,mt3,mt4],7,2);
end