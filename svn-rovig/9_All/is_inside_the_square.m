function val=is_inside_the_square(bounding_box,point)

val=(point(1)>=bounding_box(2,1) &&point(1)<=bounding_box(1,1) && point(2)>=bounding_box(2,2) &&point(2)<=bounding_box(1,2));

end