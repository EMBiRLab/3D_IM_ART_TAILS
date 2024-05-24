function out1 = five_link_2JDOF_selfcollision_con_3DOF_float_varied_len(in1,in2)
%five_link_2JDOF_selfcollision_con_3DOF_float_varied_len
%    OUT1 = five_link_2JDOF_selfcollision_con_3DOF_float_varied_len(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 9.1.
%    24-Jan-2024 15:15:26

l1 = in2(1,:);
l2 = in2(2,:);
l3 = in2(3,:);
l4 = in2(4,:);
l5 = in2(5,:);
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
q12 = in1(12,:);
q13 = in1(13,:);
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
t13 = cos(q12);
t14 = cos(q13);
t15 = sin(q1);
t16 = sin(q2);
t17 = sin(q3);
t18 = sin(q4);
t19 = sin(q5);
t20 = sin(q6);
t21 = sin(q7);
t22 = sin(q8);
t23 = sin(q9);
t24 = sin(q10);
t25 = sin(q11);
t26 = sin(q12);
t27 = sin(q13);
t28 = t2.*t4;
t29 = t2.*t17;
t30 = t4.*t15;
t31 = t5.*t16;
t32 = t15.*t17;
t34 = t16.*t18.*t19;
t35 = t2.*t3.*t5;
t36 = t3.*t4.*t6;
t38 = t3.*t5.*t15;
t39 = t3.*t4.*t19;
t42 = t3.*t17.*t18;
t43 = t6.*t16.*t18;
t48 = (t3.*t17)./2.0;
t50 = (t3.*t17)./6.0;
t52 = t3.*t17.*(5.0./6.0);
t57 = t2.*t3.*t6.*t18;
t58 = t3.*t5.*t6.*t17;
t59 = t2.*t3.*t18.*t19;
t60 = t3.*t6.*t15.*t18;
t61 = t3.*t5.*t17.*t19;
t63 = t3.*t15.*t18.*t19;
t33 = t16.*t32;
t37 = t16.*t28;
t40 = t16.*t29;
t41 = t16.*t30;
t44 = t28./2.0;
t45 = t28./6.0;
t46 = t28.*(5.0./6.0);
t47 = t30./2.0;
t49 = t30./6.0;
t51 = t30.*(5.0./6.0);
t62 = -t35;
t65 = -t43;
t72 = -t61;
t73 = -t63;
t76 = t31+t42;
t53 = -t44;
t54 = -t45;
t55 = -t46;
t56 = -t33;
t64 = -t37;
t66 = t33./2.0;
t67 = t33./6.0;
t68 = t33.*(5.0./6.0);
t69 = t40./2.0;
t70 = t40./6.0;
t71 = t40.*(5.0./6.0);
t74 = t29+t41;
t75 = t30+t40;
t80 = t7.*t76;
t85 = t8.*t20.*t76;
t86 = t20.*t21.*t76;
t96 = t39+t58+t65;
t97 = t34+t36+t72;
t77 = t28+t56;
t78 = t32+t64;
t79 = t6.*t74;
t81 = t18.*t75;
t82 = t19.*t74;
t83 = t5.*t6.*t75;
t84 = t5.*t19.*t75;
t92 = -t85;
t99 = l1.*t96;
t102 = t20.*t96;
t103 = t8.*t97;
t104 = t21.*t97;
t106 = t7.*t8.*t96;
t109 = t7.*t21.*t96;
t87 = t6.*t78;
t88 = t18.*t77;
t89 = t19.*t78;
t90 = t5.*t6.*t77;
t91 = t5.*t19.*t77;
t98 = t62+t81;
t100 = -t7.*(t35-t81);
t105 = -t8.*t20.*(t35-t81);
t108 = -t20.*t21.*(t35-t81);
t111 = t20.*t21.*(t35-t81);
t112 = -t109;
t113 = t80+t102;
t138 = t92+t104+t106;
t93 = -t89;
t94 = -t90;
t95 = t38+t88;
t114 = t59+t84+t87;
t117 = t73+t79+t91;
t118 = t9.*t113;
t123 = t10.*t22.*t113;
t124 = t22.*t23.*t113;
t139 = t86+t103+t112;
t140 = l2.*t138;
t143 = t22.*t138;
t148 = t9.*t10.*t138;
t149 = t9.*t23.*t138;
t101 = t7.*t95;
t107 = t8.*t20.*t95;
t110 = t20.*t21.*t95;
t115 = t57+t83+t93;
t116 = t60+t82+t94;
t119 = t8.*t114;
t120 = t21.*t114;
t127 = t8.*t117;
t128 = t21.*t117;
t132 = -t123;
t144 = t10.*t139;
t145 = t23.*t139;
t153 = -t149;
t154 = t118+t143;
t121 = l1.*t115;
t122 = l1.*t116;
t125 = t20.*t115;
t126 = t20.*t116;
t129 = t7.*t8.*t115;
t130 = t7.*t8.*t116;
t131 = t7.*t21.*t115;
t133 = t7.*t21.*t116;
t135 = -t127;
t158 = t11.*t154;
t162 = t12.*t24.*t154;
t163 = t24.*t25.*t154;
t180 = t132+t145+t148;
t181 = t124+t144+t153;
t134 = -t126;
t136 = -t129;
t137 = t100+t125;
t142 = t9.*(t125-t7.*(t35-t81));
t146 = t10.*t22.*(t125-t7.*(t35-t81));
t147 = t22.*t23.*(t125-t7.*(t35-t81));
t156 = t111+t119+t131;
t157 = t107+t128+t130;
t159 = -l2.*(-t120+t129+t8.*t20.*(t35-t81));
t161 = t110+t133+t135;
t164 = -t22.*(-t120+t129+t8.*t20.*(t35-t81));
t168 = l2.*(-t120+t129+t8.*t20.*(t35-t81));
t169 = -t9.*t10.*(-t120+t129+t8.*t20.*(t35-t81));
t170 = -t9.*t23.*(-t120+t129+t8.*t20.*(t35-t81));
t171 = -t162;
t173 = t22.*(-t120+t129+t8.*t20.*(t35-t81));
t183 = l3.*t180;
t184 = t24.*t180;
t185 = t12.*t181;
t186 = t25.*t181;
t189 = t11.*t12.*t180;
t190 = t11.*t25.*t180;
t141 = t101+t134;
t155 = t105+t120+t136;
t160 = l2.*t157;
t165 = t10.*t156;
t166 = t23.*t156;
t167 = t22.*t157;
t172 = t9.*t10.*t157;
t175 = t9.*t23.*t157;
t177 = t10.*t161;
t178 = t23.*t161;
t182 = t142+t173;
t195 = -t190;
t199 = t158+t184;
t221 = t171+t186+t189;
t150 = t9.*t141;
t151 = t10.*t22.*t141;
t152 = t22.*t23.*t141;
t174 = -t165;
t176 = -t167;
t179 = -t178;
t188 = t11.*t182;
t191 = t12.*t24.*t182;
t192 = t24.*t25.*t182;
t200 = t146+t166+t169;
t204 = t14.*t26.*t199;
t207 = -t12.*(-t147+t165+t9.*t23.*(-t120+t129+t8.*t20.*(t35-t81)));
t208 = -t25.*(-t147+t165+t9.*t23.*(-t120+t129+t8.*t20.*(t35-t81)));
t213 = t25.*(-t147+t165+t9.*t23.*(-t120+t129+t8.*t20.*(t35-t81)));
t222 = t163+t185+t195;
t223 = l4.*t221;
t226 = t13.*t14.*t221;
t187 = t150+t176;
t193 = -t188;
t201 = l3.*t200;
t202 = t147+t170+t174;
t203 = t24.*t200;
t205 = t152+t175+t177;
t209 = -t204;
t210 = t11.*t12.*t200;
t211 = t11.*t25.*t200;
t212 = t151+t172+t179;
t224 = t27.*t222;
t194 = t11.*t187;
t196 = t12.*t24.*t187;
t197 = t24.*t25.*t187;
t206 = -t201;
t214 = t12.*t205;
t215 = t25.*t205;
t216 = l3.*t212;
t217 = t24.*t212;
t219 = t11.*t12.*t212;
t220 = t11.*t25.*t212;
t225 = t193+t203;
t227 = -t14.*t26.*(t188-t203);
t230 = t192+t207+t211;
t231 = t191+t210+t213;
t243 = t209+t224+t226;
t198 = -t194;
t218 = -t215;
t229 = -t14.*t26.*(t194-t217);
t232 = t27.*t230;
t233 = l4.*t231;
t235 = t13.*t14.*t231;
t237 = t197+t214+t220;
t244 = l5.*t243;
t228 = t198+t217;
t234 = -t233;
t236 = -t235;
t238 = t196+t218+t219;
t239 = t27.*t237;
t246 = -l5.*(-t232+t235+t14.*t26.*(t188-t203));
t240 = l4.*t238;
t241 = t13.*t14.*t238;
t245 = t227+t232+t236;
t242 = -t241;
t248 = -l5.*(-t239+t241+t14.*t26.*(t194-t217));
t249 = l5.*(-t239+t241+t14.*t26.*(t194-t217));
t247 = t229+t239+t242;
mt1 = [-(t51+t71+t121+t168+t206).^2-(t55+t68+t122+t160+t216).^2-(t52+t99+t140+t183).^2+9.388888888888889e-2;-(t47+t69+t121+t168+t206).^2-(t53+t66+t122+t160+t216).^2-(t48+t99+t140+t183).^2+9.388888888888889e-2;-(t49+t70+t121+t168+t206).^2-(t54+t67+t122+t160+t216).^2-(t50+t99+t140+t183).^2+9.388888888888889e-2;-(t52+t99+t140+t183+t223).^2-(t51+t71+t121+t168+t206+t234).^2-(t55+t68+t122+t160+t216+t240).^2+9.388888888888889e-2];
mt2 = [-(t48+t99+t140+t183+t223).^2-(t47+t69+t121+t168+t206+t234).^2-(t53+t66+t122+t160+t216+t240).^2+9.388888888888889e-2;-(t50+t99+t140+t183+t223).^2-(t49+t70+t121+t168+t206+t234).^2-(t54+t67+t122+t160+t216+t240).^2+9.388888888888889e-2;-(t52+t99+t140+t183+t223+t244).^2-(t51+t71+t121+t168+t206+t234+t246).^2-(t55+t68+t122+t160+t216+t240+t249).^2+9.388888888888889e-2;-(t48+t99+t140+t183+t223+t244).^2-(t47+t69+t121+t168+t206+t234+t246).^2-(t53+t66+t122+t160+t216+t240+t249).^2+9.388888888888889e-2];
mt3 = [-(t50+t99+t140+t183+t223+t244).^2-(t49+t70+t121+t168+t206+t234+t246).^2-(t54+t67+t122+t160+t216+t240+t249).^2+9.388888888888889e-2];
out1 = [mt1;mt2;mt3];
