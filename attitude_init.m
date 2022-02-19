% This function initializes attitudes for a MPF-initialization.

function attitudes = attitude_init(n1,n2,n3,smph_roll,smph_pitch)

attitudes = zeros(3,n1*n2*n3);

s1 = 0.2; % Variation of smartphone roll.
s2 = 0.2; % Variation of smartphone pitch.

if(n1~=0)
    roll_ser = (-s1:2*s1/(n1-1):s1)+smph_roll; 
else
    roll_ser = smph_roll;
end

if(n2~=0)
    pitch_ser = (-s2:2*s2/(n2-1):s2)+smph_pitch; 
else
    pitch_ser = smph_pitch;
end

[mat1,mat2,mat3] = meshgrid(roll_ser,pitch_ser,2*pi/n3:2*pi/n3:2*pi);

for particle = 1:n1*n2*n3
    attitudes(1,particle) = mat1(particle);
    attitudes(2,particle) = mat2(particle);
    attitudes(3,particle) = mat3(particle);
end
    
end

