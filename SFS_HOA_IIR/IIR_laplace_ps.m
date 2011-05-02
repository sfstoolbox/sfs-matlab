function [drivFuncCoeff] = IIR_laplace_ps(config)

% Parameters
R=1.5; % Radius of LS-Array
c=343; % speed of sound
f=linspace(1,20000,300); % frequency vector

% Position of sound source (relative to center of LS-array)
xs_rel = ( config.xps(1) - config.xref(1) );
ys_rel = ( config.xps(2) - config.xref(2) );

rs = sqrt(  xs_rel.^2 + ys_rel.^2 ); %  get distance from point source
alphas = atan2( ys_rel, xs_rel );   % get angle from position of point source


% get alpha_0 (from config declared in HOA_main.m)
alpha_0=atan2( config.x0( 2,: ), config.x0( 1,: ) );



% allocate memory for driving functions
drivFuncCoeff = zeros(56,length(f));

% allocate memory for inverse laplace transformed filter
invLaplace = zeros(56,length(f));


for j = 1:56 % loop over all loudpspeakers
    disp([j 56]) %display calculation (time) for every LS
    
    i = 1; % set fill-index of invLaplace
    
    for m = -27:28 % band limitation for loudspeakers
    
        B=zeros(1,abs(m)+1); % allocate memory for denominator of driving function
        A=zeros(1,abs(m)+1); % allocate memory for numerator of driving function

    
        for k=0:abs(m) % calculate inner sum of denominator and numerator
           
            B(k+1) = beta1(abs(m),k);
            A(k+1) = B(k+1).*(rs./c)^k.*(R/rs)^(abs(m)+1);   
            B(k+1) = B(k+1).* (R./c)^k;
           
        end
        
        
    B = fliplr(B); % flip vector entries for descending order of denominator
    A = fliplr(A);
    
    % calculate frequency response for single filter
    invLaplace(i,:)=freqs(A,B,2*pi*f);

    % multiply with remaining terms 
    invLaplace(i,:)= invLaplace(i,:).*(1./(2*pi*R)) .* (exp(+1i*2*pi*f/c.*(R-rs))) .*...
                                           exp(1i.*m.* ( alpha_0(j) - alphas ) );
    
     i=i+1; % set index up to fill next line of invLaplace

    end

% get the number of columns in invLaplace (number of loudspeakers) 
u = size(invLaplace);
p = u(1);
    
        for q = 1:p % add all calculated filters to get driving functions
    
            drivFuncCoeff(j,:) = invLaplace(q,:) + drivFuncCoeff(j,:); 
    
    
        end
         
    

end



end

