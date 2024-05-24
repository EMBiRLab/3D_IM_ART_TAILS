function out1 = four_link_2JDOF_selfcollision_con_3DOF_float(in1)
%four_link_2JDOF_selfcollision_con_3DOF_float
%    OUT1 = four_link_2JDOF_selfcollision_con_3DOF_float(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    21-Jan-2024 14:50:20

q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
q4 = in1(4,:);
q5 = in1(5,:);
q6 = in1(6,:);
q7 = in1(7,:);
q8 = in1(8,:);
q9 = in1(9,:);
q10 = in1(10,:);
q11 = in1(11,:);
t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = cos(q4);
t6 = cos(q5);
t7 = cos(q6);
t8 = cos(q7);
t9 = cos(q8);
t10 = cos(q9);
t11 = cos(q10);
t12 = cos(q11);
t13 = sin(q1);
t14 = sin(q2);
t15 = sin(q3);
t16 = sin(q4);
t17 = sin(q5);
t18 = sin(q6);
t19 = sin(q7);
t20 = sin(q8);
t21 = sin(q9);
t22 = sin(q10);
t23 = sin(q11);
t24 = t2.*t4;
t25 = t2.*t15;
t26 = t4.*t13;
t27 = t5.*t14;
t28 = t13.*t15;
t30 = t14.*t16.*t17;
t31 = t2.*t3.*t5;
t32 = t3.*t4.*t6;
t34 = t3.*t5.*t13;
t35 = t3.*t4.*t17;
t38 = t3.*t15.*t16;
t39 = t6.*t14.*t16;
t44 = (t3.*t15)./2.0;
t46 = (t3.*t15)./6.0;
t48 = t3.*t15.*(5.0./6.0);
t56 = t2.*t3.*t6.*t16;
t57 = t3.*t5.*t6.*t15;
t58 = t2.*t3.*t16.*t17;
t59 = t3.*t6.*t13.*t16;
t60 = t3.*t5.*t15.*t17;
t62 = t3.*t13.*t16.*t17;
t29 = t14.*t28;
t33 = t14.*t24;
t36 = t14.*t25;
t37 = t14.*t26;
t40 = t24./2.0;
t41 = t24./6.0;
t42 = t24.*(5.0./6.0);
t43 = t26./2.0;
t45 = t26./6.0;
t47 = t26.*(5.0./6.0);
t61 = -t31;
t64 = -t39;
t71 = -t60;
t72 = -t62;
t73 = t35.*(3.0./8.0);
t77 = t39.*(3.0./8.0);
t81 = t27+t38;
t82 = t56.*(3.0./8.0);
t83 = t57.*(3.0./8.0);
t84 = t59.*(3.0./8.0);
t49 = -t40;
t50 = -t41;
t51 = -t42;
t52 = -t29;
t53 = -t43;
t54 = -t45;
t55 = -t47;
t63 = -t33;
t65 = t29./2.0;
t66 = t29./6.0;
t67 = t29.*(5.0./6.0);
t68 = t36./2.0;
t69 = t36./6.0;
t70 = t36.*(5.0./6.0);
t78 = -t77;
t79 = t25+t37;
t80 = t26+t36;
t85 = -t82;
t89 = t7.*t81;
t94 = t8.*t18.*t81;
t95 = t18.*t19.*t81;
t113 = t35+t57+t64;
t114 = t30+t32+t71;
t74 = -t68;
t75 = -t69;
t76 = -t70;
t86 = t24+t52;
t87 = t28+t63;
t88 = t6.*t79;
t90 = t16.*t80;
t91 = t17.*t79;
t92 = t5.*t6.*t80;
t93 = t5.*t17.*t80;
t101 = -t94;
t106 = t94.*(3.0./8.0);
t118 = t18.*t113;
t119 = t8.*t114;
t120 = t19.*t114;
t122 = t7.*t8.*t113;
t125 = t7.*t19.*t113;
t96 = t6.*t87;
t97 = t16.*t86;
t98 = t17.*t87;
t99 = t5.*t6.*t86;
t100 = t5.*t17.*t86;
t103 = t91.*(3.0./8.0);
t105 = t92.*(3.0./8.0);
t109 = -t106;
t115 = t61+t90;
t116 = -t7.*(t31-t90);
t121 = -t8.*t18.*(t31-t90);
t124 = -t18.*t19.*(t31-t90);
t127 = t18.*t19.*(t31-t90);
t128 = -t125;
t129 = t120.*(3.0./8.0);
t130 = t8.*t18.*(t31-t90).*(-3.0./8.0);
t131 = t122.*(3.0./8.0);
t133 = t89+t118;
t163 = t101+t120+t122;
t102 = -t98;
t104 = -t99;
t107 = t98.*(3.0./8.0);
t108 = -t105;
t110 = t99.*(3.0./8.0);
t112 = t34+t97;
t134 = t58+t93+t96;
t137 = t72+t88+t100;
t138 = t9.*t133;
t141 = t10.*t20.*t133;
t142 = t20.*t21.*t133;
t164 = t95+t119+t128;
t167 = t20.*t163;
t172 = t9.*t10.*t163;
t173 = t9.*t21.*t163;
t111 = -t110;
t117 = t7.*t112;
t123 = t8.*t18.*t112;
t126 = t18.*t19.*t112;
t135 = t56+t92+t102;
t136 = t59+t91+t104;
t139 = t8.*t134;
t140 = t19.*t134;
t145 = t8.*t137;
t146 = t19.*t137;
t150 = -t141;
t156 = t141.*(3.0./8.0);
t168 = t10.*t164;
t169 = t21.*t164;
t177 = -t173;
t180 = t172.*(3.0./8.0);
t182 = t138+t167;
t132 = t123.*(3.0./8.0);
t143 = t18.*t135;
t144 = t18.*t136;
t147 = t7.*t8.*t135;
t148 = t7.*t8.*t136;
t149 = t7.*t19.*t135;
t151 = t7.*t19.*t136;
t153 = -t145;
t154 = t140.*(3.0./8.0);
t157 = t146.*(3.0./8.0);
t158 = -t156;
t178 = t169.*(3.0./8.0);
t200 = t12.*t22.*t182.*(3.0./8.0);
t208 = t150+t169+t172;
t209 = t142+t168+t177;
t152 = -t144;
t155 = -t147;
t159 = t147.*(3.0./8.0);
t160 = t148.*(3.0./8.0);
t162 = t116+t143;
t166 = t9.*(t143-t7.*(t31-t90));
t170 = t10.*t20.*(t143-t7.*(t31-t90));
t171 = t20.*t21.*(t143-t7.*(t31-t90));
t184 = t127+t139+t149;
t185 = t123+t146+t148;
t186 = t126+t151+t153;
t187 = -t20.*(-t140+t147+t8.*t18.*(t31-t90));
t191 = -t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90));
t192 = -t9.*t21.*(-t140+t147+t8.*t18.*(t31-t90));
t194 = t20.*(-t140+t147+t8.*t18.*(t31-t90));
t203 = -t200;
t204 = t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90)).*(-3.0./8.0);
t212 = t23.*t209.*(3.0./8.0);
t213 = t11.*t12.*t208.*(3.0./8.0);
t161 = -t159;
t165 = t117+t152;
t179 = t170.*(3.0./8.0);
t183 = t121+t140+t155;
t188 = t10.*t184;
t189 = t21.*t184;
t190 = t20.*t185;
t193 = t9.*t10.*t185;
t196 = t9.*t21.*t185;
t198 = t10.*t186;
t199 = t21.*t186;
t210 = t166+t194;
t174 = t9.*t165;
t175 = t10.*t20.*t165;
t176 = t20.*t21.*t165;
t195 = -t188;
t197 = -t190;
t201 = t189.*(3.0./8.0);
t202 = -t199;
t205 = t193.*(3.0./8.0);
t206 = t199.*(3.0./8.0);
t214 = t12.*t22.*t210.*(3.0./8.0);
t216 = t170+t189+t191;
t220 = t23.*(-t171+t188+t9.*t21.*(-t140+t147+t8.*t18.*(t31-t90))).*(-3.0./8.0);
t222 = t23.*(-t171+t188+t9.*t21.*(-t140+t147+t8.*t18.*(t31-t90))).*(3.0./8.0);
t181 = t175.*(3.0./8.0);
t207 = -t206;
t211 = t174+t197;
t217 = t171+t192+t195;
t218 = t176+t196+t198;
t219 = t175+t193+t202;
t221 = t11.*t12.*t216.*(3.0./8.0);
t215 = t12.*t22.*t211.*(3.0./8.0);
t223 = t23.*t218.*(3.0./8.0);
t225 = t11.*t12.*t219.*(3.0./8.0);
t224 = -t223;
mt1 = [-(t47+t70+t82+t105-t107-t154+t159-t179-t201+t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90)).*(3.0./8.0)+t8.*t18.*(t31-t90).*(3.0./8.0)).^2-(t48+t73+t78+t83+t109+t129+t131+t158+t178+t180).^2-(t51+t67+t84+t103+t111+t132+t157+t160+t181+t205+t207).^2+9.388888888888889e-2];
mt2 = [-(t43+t68+t82+t105-t107-t154+t159-t179-t201+t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90)).*(3.0./8.0)+t8.*t18.*(t31-t90).*(3.0./8.0)).^2-(t44+t73+t78+t83+t109+t129+t131+t158+t178+t180).^2-(t49+t65+t84+t103+t111+t132+t157+t160+t181+t205+t207).^2+9.388888888888889e-2];
mt3 = [-(t45+t69+t82+t105-t107-t154+t159-t179-t201+t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90)).*(3.0./8.0)+t8.*t18.*(t31-t90).*(3.0./8.0)).^2-(t46+t73+t78+t83+t109+t129+t131+t158+t178+t180).^2-(t50+t66+t84+t103+t111+t132+t157+t160+t181+t205+t207).^2+9.388888888888889e-2];
mt4 = [-(t51+t67+t84+t103+t111+t132+t157+t160+t181+t205+t207+t215+t224+t225).^2-(t47+t70+t82+t105-t107-t154+t159-t179-t201-t214+t220-t221+t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90)).*(3.0./8.0)+t8.*t18.*(t31-t90).*(3.0./8.0)).^2-(t48+t73+t78+t83+t109+t129+t131+t158+t178+t180+t203+t212+t213).^2+9.388888888888889e-2];
mt5 = [-(t49+t65+t84+t103+t111+t132+t157+t160+t181+t205+t207+t215+t224+t225).^2-(t43+t68+t82+t105-t107-t154+t159-t179-t201-t214+t220-t221+t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90)).*(3.0./8.0)+t8.*t18.*(t31-t90).*(3.0./8.0)).^2-(t44+t73+t78+t83+t109+t129+t131+t158+t178+t180+t203+t212+t213).^2+9.388888888888889e-2];
mt6 = [-(t50+t66+t84+t103+t111+t132+t157+t160+t181+t205+t207+t215+t224+t225).^2-(t45+t69+t82+t105-t107-t154+t159-t179-t201-t214+t220-t221+t9.*t10.*(-t140+t147+t8.*t18.*(t31-t90)).*(3.0./8.0)+t8.*t18.*(t31-t90).*(3.0./8.0)).^2-(t46+t73+t78+t83+t109+t129+t131+t158+t178+t180+t203+t212+t213).^2+9.388888888888889e-2];
out1 = [mt1;mt2;mt3;mt4;mt5;mt6];
