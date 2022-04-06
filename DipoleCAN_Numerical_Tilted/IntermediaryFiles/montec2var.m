function I=montec2var(g,a,b,c,d,n)


%{
%I=montec2var('g',a,b,c,d,n)
%Calculates doble integral using Montecarlo method (fixed boundaries)
I = zeros( n, 3 );
for k=1:n
    s=2.^k;
    t=rand(2,s);
    x=a+t(1,:).*(b-a);			%extremos a,b fijos
    y=c+t(2,:).*(d-c);			%extremos c,d fijos
    
    z = zeros( s, 3 );
    for i=1:s
        z(i, :)=feval(g,x(i),y(i));
    end

    I(k,1)=(((b-a).*(d-c))./s).*sum(z(:,1));
    I(k,2)=(((b-a).*(d-c))./s).*sum(z(:,2));
    I(k,3)=(((b-a).*(d-c))./s).*sum(z(:,3));
end

hold on
    plot(I(:, 1));
    plot(I(:, 2));
    plot(I(:, 3));
%}
    
    
    
    
    t1=rand(1,n);
    t2=rand(1,n);
    x=a+t1(1,:).*(b-a);
    y=c+t2(1,:).*(d-c);
    p = 1./((b-a)*(d-c));
    
    Sum = zeros( 1, 3 );
    z = zeros( n, 3 );
    for i=1:n
        z(i, :) = feval(g,x(i),y(i));
    end
        Sum(:,1) = sum(z(:, 1));
        Sum(:,2) = sum(z(:, 2));
        Sum(:,3) = sum(z(:, 3));
        
        I = Sum ./ (p*n);
 
end
    