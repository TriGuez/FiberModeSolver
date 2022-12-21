function [neff, LP] = ModeSolver(RImap, x, y, varargin)
    warning('off','all');
    nModes = 10;
    coreRadius = 5e-6;
    lambda = 1064e-9;
    plot = true ;
    target = 0;
    IndexContour = true;
    for ii = 1:2:numel(varargin)
        switch(lower(varargin{ii}))
            case 'nmodes'
                nModes = varargin{ii+1};
            case 'coreradius'
                coreRadius = varargin{ii+1};
            case 'lambda'
                lambda = varargin{ii+1};
            case 'plot'
                plot = varargin{ii+1};
            case 'target'
                target = varargin{ii+1};
            case 'indexcontour'
                IndexContour = varargin{ii+1};
            otherwise
                error('Unknown argument ''%s'' ', varargin{ii})
        end
    end
    N = length(x);
    h = max(x(2)-x(1), y(2)-y(1));
    k0 = 2*pi/lambda;
    target = (target*k0).^2;
    lowerdiag = ones(N^2,1);
    lowerdiag(N:N:end) = 0;
    upperdiag = circshift(lowerdiag,1);
    Op = ( spdiags(-4/h.^2.*ones(N^2,1)+k0.^2.*((RImap(:)).^2),0,N^2,N^2)+...
        spdiags((1/h^2).*lowerdiag,-1,N^2,N^2) + spdiags((1/h^2).*upperdiag,1,N^2,N^2) + ...
        spdiags((1/h^2).*ones(N^2,1),-N,N^2,N^2) + spdiags((1/h^2).*ones(N^2,1),N,N^2,N^2));
    if target ~= 0
        [V,D] = eigs(Op, nModes,target);
    else
        [V,D] = eigs(Op, nModes, 'largestreal');
    end
    D = sqrt(diag(D)/k0.^2);
    idxs = find(~isnan(D));
    D = D(idxs);
    V = reshape(V(:,idxs), N,N,length(idxs));
    I = abs(V).^2;
    if isempty(D)
        error('Solver did not converge, 0 mode found')
    else
        fprintf(1,'%d converged mode(s) found, processing core modes\n', length(D))
        overlaps = zeros(1,length(idxs));
        idx_x_core = find(abs(x) < coreRadius);
        idx_y_core = find(abs(y) < coreRadius);
        for jk = 1:length(D)
            overlaps(jk) = trapz(y(idx_y_core), trapz(x(idx_x_core),I(idx_x_core,idx_y_core,jk)))./trapz(y,trapz(x,(I(:,:,jk))));
        end
        coreModeIdxs = find(overlaps > 0.7);
        if isempty(coreModeIdxs)
            fprintf(1,'No core mode found\n')
        else
            fprintf(1,'%d core mode(s) found\n', length(coreModeIdxs))
        end
        LP = V(:,:,coreModeIdxs);
        neff = D(coreModeIdxs);
    end
    if plot
        x1 = linspace(min(x),max(x),2048);
        y1 = linspace(min(y),max(y),2048);
        [X,Y] = meshgrid(x,y);
        [X1,Y1] = meshgrid(x1,y1);
        for jk = 1:length(neff)
            Iplot = abs(LP(:,:,jk)).^2;
            Iplot = Iplot./(max(max(Iplot)));
            Iplot = interp2(X,Y,Iplot,X1,Y1,'spline');
            figure(100+jk)
            imagesc(x.*1e6,y.*1e6,Iplot)
            colormap jet
            if IndexContour
                cRange = caxis;
                hold on
                contour(x.*1e6,y.*1e6,real(RImap),'-white')
                caxis(cRange);
            end
            xlabel('x (\mu m)')
            ylabel('y (\mu m)')
            title(['n_{eff} = ' num2str(neff(jk),16)])
        end
    end
end
    
    