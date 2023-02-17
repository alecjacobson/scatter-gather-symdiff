function [sym_f,sym_dfdx,sym_d2fdx2] = spring_3d_sym(in1,c_11,c_21)
%SPRING_3D_SYM
%    [SYM_F,SYM_DFDX,SYM_D2FDX2] = SPRING_3D_SYM(IN1,C_11,C_21)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    14-Jul-2022 20:53:57

x1 = in1(:,1);
x2 = in1(:,2);
x3 = in1(:,3);
x4 = in1(:,4);
x5 = in1(:,5);
x6 = in1(:,6);
t2 = x1.*2.0;
t3 = x2.*2.0;
t4 = x3.*2.0;
t5 = x4.*2.0;
t6 = x5.*2.0;
t7 = x6.*2.0;
t8 = -x2;
t10 = -x4;
t12 = -x6;
t9 = -t3;
t11 = -t5;
t13 = -t7;
t14 = t8+x1;
t15 = t10+x3;
t16 = t12+x5;
t17 = t2+t9;
t18 = t4+t11;
t19 = t6+t13;
t20 = t14.^2;
t21 = t15.^2;
t22 = t16.^2;
t23 = t17.^2;
t24 = t18.^2;
t25 = t19.^2;
t26 = t20+t21+t22;
t27 = 1.0./t26;
t28 = sqrt(t26);
t29 = 1.0./t28;
t31 = -t28;
t33 = (c_11.*t23.*t27)./2.0;
t34 = (c_11.*t24.*t27)./2.0;
t35 = (c_11.*t25.*t27)./2.0;
t39 = (c_11.*t17.*t18.*t27)./2.0;
t40 = (c_11.*t17.*t19.*t27)./2.0;
t41 = (c_11.*t18.*t19.*t27)./2.0;
t30 = t29.^3;
t32 = c_21+t31;
sym_f = c_11.*t32.^2;
if nargout > 1
    t36 = -t33;
    t37 = -t34;
    t38 = -t35;
    t42 = -t39;
    t43 = -t40;
    t44 = -t41;
    t45 = c_11.*t29.*t32.*2.0;
    t47 = c_11.*t17.*t29.*t32;
    t48 = c_11.*t18.*t29.*t32;
    t49 = c_11.*t19.*t29.*t32;
    sym_dfdx = [-t47;t47;-t48;t48;-t49;t49];
end
if nargout > 2
    t50 = (c_11.*t23.*t30.*t32)./2.0;
    t51 = (c_11.*t24.*t30.*t32)./2.0;
    t52 = (c_11.*t25.*t30.*t32)./2.0;
    t56 = (c_11.*t17.*t18.*t30.*t32)./2.0;
    t57 = (c_11.*t17.*t19.*t30.*t32)./2.0;
    t58 = (c_11.*t18.*t19.*t30.*t32)./2.0;
    t46 = -t45;
    t53 = -t50;
    t54 = -t51;
    t55 = -t52;
    t59 = -t56;
    t60 = -t57;
    t61 = -t58;
    t62 = t39+t56;
    t63 = t40+t57;
    t64 = t41+t58;
    t65 = t42+t59;
    t66 = t43+t60;
    t67 = t44+t61;
    t68 = t33+t46+t50;
    t69 = t34+t46+t51;
    t70 = t35+t46+t52;
    t71 = t36+t45+t53;
    t72 = t37+t45+t54;
    t73 = t38+t45+t55;
    sym_d2fdx2 = [t68;t71;t62;t65;t63;t66;t71;t68;t65;t62;t66;t63;t62;t65;t69;t72;t64;t67;t65;t62;t72;t69;t67;t64;t63;t66;t64;t67;t70;t73;t66;t63;t67;t64;t73;t70];
end