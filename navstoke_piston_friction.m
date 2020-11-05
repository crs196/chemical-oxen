syms x h(x,t) p(x,t) u mu t
lhs=diff(h^3*diff(p,x),x);
rhs=6*u*mu*diff(h,x)+12*mu*diff(h,t);
dsolve(lhs==rhs);