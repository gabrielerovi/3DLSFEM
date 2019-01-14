function A=area(node) 

   x1=node(1,1);x2=node(2,1);x3=node(3,1);
   y1=node(1,2);y2=node(2,2);y3=node(3,2);
   
    A=abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2.0;

end