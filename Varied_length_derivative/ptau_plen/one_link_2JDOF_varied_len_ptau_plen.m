function ptau_plen = one_link_2JDOF_varied_len_ptau_plen(in1,in2,in3,l1)
%one_link_2JDOF_varied_len_ptau_plen
%    PTAU_PLEN = one_link_2JDOF_varied_len_ptau_plen(IN1,IN2,IN3,L1)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    24-Jan-2024 11:17:42

q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
qd1 = in2(1,:);
qd2 = in2(2,:);
qd3 = in2(3,:);
qd4 = in2(4,:);
qd5 = in2(5,:);
qdd1 = in3(1,:);
qdd2 = in3(2,:);
qdd3 = in3(3,:);
qdd4 = in3(4,:);
qdd5 = in3(5,:);
t2 = cos(q2);
t3 = cos(q3);
t4 = cos(q4);
t5 = cos(q5);
t6 = sin(q2);
t7 = sin(q3);
t8 = sin(q4);
t9 = sin(q5);
t10 = l1.^2;
t21 = qd3./2.0;
t22 = qdd3./2.0;
t11 = qd2.*t3;
t12 = qdd1.*t2;
t13 = qdd2.*t3;
t14 = qd1.*t6;
t15 = qd2.*t7;
t16 = qdd1.*t6;
t17 = qdd2.*t7;
t18 = qd1.*qd2.*t2;
t20 = t10.*2.0;
t23 = qd1.*t2.*t3;
t25 = qd1.*t2.*t7;
t19 = qd2.*t14;
t24 = qd3+t14;
t26 = -t13;
t28 = t14./2.0;
t29 = t16./2.0;
t30 = t18./2.0;
t31 = -t25;
t34 = t15+t23;
t35 = qdd3+t16+t18;
t37 = t20+1.0./6.0e+2;
t27 = -t19;
t32 = t4.*t24;
t33 = t8.*t24;
t38 = qd4+t34;
t39 = qd3.*t34;
t40 = t11+t31;
t41 = t21+t28;
t42 = t4.*t34;
t43 = t4.*t35;
t44 = t8.*t35;
t47 = qd4.*t8.*t34;
t54 = t8.*t9.*t34;
t65 = (t5.*t8.*t34)./2.0;
t68 = t22+t29+t30;
t36 = t12+t27;
t45 = l1.*t42;
t46 = t5.*t38;
t48 = t9.*t38;
t49 = qd3.*t40;
t52 = t5.*t41;
t53 = t9.*t41;
t55 = t4.*t40;
t56 = t8.*t40;
t58 = -t47;
t60 = (qd4.*t42)./2.0;
t61 = t47./2.0;
t62 = -t54;
t66 = t54./2.0;
t69 = t5.*t68;
t71 = t9.*t68.*2.0;
t50 = t3.*t36;
t51 = t7.*t36;
t57 = t52.*2.0;
t59 = -t48;
t63 = -t56;
t64 = -t61;
t67 = -t66;
t70 = t69.*2.0;
t72 = -t69;
t74 = -t71;
t75 = t33+t55;
t88 = t53+t65;
t73 = -t70;
t76 = qd4.*t75;
t77 = t32+t63;
t78 = t5.*t75;
t79 = t9.*t75;
t86 = t26+t39+t51;
t87 = t17+t49+t50;
t90 = t52+t67;
t91 = qd5.*t88;
t80 = qd5+t77;
t81 = qd4.*t77;
t82 = -t76;
t89 = qdd4+t87;
t92 = t4.*t86;
t93 = t4.*t87;
t94 = t8.*t86;
t95 = t91.*2.0;
t97 = qd5.*t90.*2.0;
t99 = l1.*t90.*2.0;
t102 = (t8.*t87)./2.0;
t105 = t46+t79;
t107 = t59+t78;
t109 = -qd5.*(t48-t78);
t115 = -t42.*(t48-t78);
t117 = -t45.*(t48-t78);
t141 = t90.*(t48-t78).*-2.0;
t83 = t80.^2;
t84 = l1.*t80.*2.0;
t96 = t5.*t89;
t98 = -t92;
t100 = -t97;
t101 = t93./2.0;
t103 = t37.*t80;
t104 = (t9.*t89)./3.0e+2;
t106 = t105.^2;
t108 = l1.*t105.*2.0;
t111 = (qd5.*t105)./3.0e+2;
t112 = t37.*t105;
t113 = t42.*t105;
t114 = t60+t102;
t124 = t80.*t88.*2.0;
t127 = t80.*t90.*2.0;
t129 = l1.*t80.*t88.*-2.0;
t131 = l1.*t80.*(t48-t78).*-2.0;
t132 = t80.*(t48-t78).*(-1.0./3.0e+2);
t133 = t88.*t105.*2.0;
t134 = (t80.*(t48-t78))./3.0e+2;
t136 = qdd5+t43+t82+t94;
t143 = l1.*t141;
t149 = l1.*t105.*(t48-t78).*-2.0;
t151 = t105.*(t48-t78).*(-1.0./3.0e+2);
t153 = (t105.*(t48-t78))./3.0e+2;
t85 = l1.*t83.*2.0;
t110 = l1.*t106.*2.0;
t116 = t64+t101;
t118 = t42+t108;
t119 = t9.*t114;
t120 = l1.*(t61-t101).*-2.0;
t121 = t5.*t114.*2.0;
t125 = t57+t62+t84;
t126 = t84.*t88;
t128 = -t124;
t130 = t45+t112;
t135 = t99+t103;
t137 = t88.*t108;
t138 = t44+t81+t98;
t142 = l1.*t136.*2.0;
t146 = -t103.*(t48-t78);
t150 = t37.*t136;
t156 = -t112.*(t48-t78);
t122 = t119.*2.0;
t123 = -t121;
t139 = (t42.*t125)./2.0;
t144 = -t142;
t145 = t9.*t138;
t147 = (t5.*t138)./3.0e+2;
t152 = -t150;
t154 = t90.*t118;
t155 = t80.*t130;
t157 = t72+t91+t119;
t158 = t105.*t135;
t140 = -t139;
t148 = -t147;
t159 = l1.*t157.*2.0;
t160 = -t158;
t161 = t96+t109+t145;
t164 = t74+t85+t100+t110+t113+t123+t127;
t167 = t73+t95+t115+t122+t128+t144+t149;
t168 = -t9.*(t70-t95-t122+t124+t142+t42.*(t48-t78)+t108.*(t48-t78));
t169 = t9.*(t70-t95-t122+t124+t142+t42.*(t48-t78)+t108.*(t48-t78));
t170 = t5.*(t70-t95-t122+t124+t142+t42.*(t48-t78)+t108.*(t48-t78)).*(-1.0./2.0);
t162 = l1.*t161.*2.0;
t163 = t37.*t161;
t165 = t5.*t164;
t166 = (t9.*t164)./2.0;
t173 = t117+t129+t152+t153+t156+t159;
t174 = -t4.*(t126+t150+t151-t159+t45.*(t48-t78)+t112.*(t48-t78));
t175 = -t8.*(t126+t150+t151-t159+t45.*(t48-t78)+t112.*(t48-t78));
t179 = t104+t111+t140+t148+t154+t155+t160;
t171 = t58+t93+t131+t133+t141+t162;
t172 = t4.*(t47-t93-t133-t162+t84.*(t48-t78)+t90.*(t48-t78).*2.0).*(-1.0./2.0);
t176 = t120+t134+t137+t143+t146+t163;
t177 = -t5.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0);
t178 = -t9.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0);
t180 = t9.*t179;
t181 = t5.*t179;
t183 = t165+t169;
t182 = -t181;
t184 = (t8.*t183)./2.0;
t187 = -t4.*(t181+t9.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0));
t188 = -t8.*(t181+t9.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0));
t185 = -t184;
t186 = t178+t182;
t189 = t175+t187;
et1 = -t2.*(t3.*(-t180+t184+(t4.*(t47-t93-t133-t162+t84.*(t48-t78)+t90.*(t48-t78).*2.0))./2.0+t5.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0))-t7.*(t4.*(t181+t9.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0))+t8.*(t126+t150+t151-t159+t45.*(t48-t78)+t112.*(t48-t78))));
et2 = -t6.*(t166+t170+t174+t8.*(t181+t9.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0)));
mt1 = [et1+et2;-t7.*(-t180+t184+(t4.*(t47-t93-t133-t162+t84.*(t48-t78)+t90.*(t48-t78).*2.0))./2.0+t5.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0))-t3.*(t4.*(t181+t9.*(t132-t163+l1.*(t61-t101).*2.0+t99.*(t48-t78)+t103.*(t48-t78)-l1.*t88.*t105.*2.0))+t8.*(t126+t150+t151-t159+t45.*(t48-t78)+t112.*(t48-t78)))];
mt2 = [-t166+t188+(t5.*(t70-t95-t122+t124+t142+t42.*(t48-t78)+t108.*(t48-t78)))./2.0+t4.*(t126+t150+t151-t159+t45.*(t48-t78)+t112.*(t48-t78));t177+t180;t126+t150+t151-t159+t45.*(t48-t78)+t112.*(t48-t78)];
ptau_plen = [mt1;mt2];
