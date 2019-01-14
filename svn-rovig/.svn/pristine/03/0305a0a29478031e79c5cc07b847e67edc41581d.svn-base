function [xn,T]=OrthogonalTransformation(x,normal)


theta=(2*pi)*normal'* x / (norm(normal)*norm(x));

T=[ cos(theta) -sin(theta);
    sin(theta)  cos(theta);];

xn=T*x;

end