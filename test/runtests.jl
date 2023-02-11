using Test
using AtiyahBott

# first test
println( "Testing O1" );
P = O1()^2;
res = AtiyahBottFormula(2,2,6-1,P);
if res[1] == 1//1
    @test true
else
    @test false
end

#second test
println( "Testing Incidency and Contact" )
d = 4;
P = Vector{Any}(undef, d+1);
for a in 0:d
    P[a+1] = Incidency(3)^a*Incidency(2)^(2*(d-a)+1)*Contact();
end
res = AtiyahBottFormula(3,d,0,P);

if res == [ 1089024, 96512, 9408, 1024, 128 ] 
    @test true
else
    @test false
end

#third test
println( "Testing Hypersurface" );
d = 4 #for other values of d, change this line
P = Hypersurface(5);
res = AtiyahBottFormula(4,d,0,P);

if res == [ 15517926796875//64 ]
    @test true 
else  
    @test false
end

#fourth test
println( "Testing R1" );
d = 4 #for other values of d, change this line
P = R1(1)^2;
res = AtiyahBottFormula(1,d,0,P);

if res == [ 1//64 ]
    @test true
else
    @test false
end