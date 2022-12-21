function MFD = ModeFieldDiameter(x,y,Field,varargin)
    
type = '4sigma';
    for ii=1:numel(varargin)
        switch varargin{ii}
            case '1/e'
                type = '1/e';
            case '1/e2'
                type = '1/e2';
            case 'FWHM'
                type = 'FWHM';
            case '4sigma'
                type = '4sigma';
            otherwise
                error('Unknown argument ''%f'' ', varargin{ii});
        end
    end
    if ~isreal(Field)
        Field = abs(Field).^2;
    end
    Field = Field./(max(max(Field)));
    FT_x = Field(:,length(Field)/2);
    FT_y = Field(length(Field)/2,:);
    switch type
        case '1/e'
            idx_x = find(FT_x > max(FT_x/(exp(1))));
            idx_y = find(FT_y > max(FT_y/(exp(1))));
            ax = min(idx_x); bx = max(idx_x);
            ay = min(idx_x); by = max(idx_y);
            MFD = max((x(bx)-x(ax)),(x(by)-x(ay)));
        case '1/e2'
            idx_x = find(FT_x > max(FT_x/(exp(2))));
            idx_y = find(FT_y > max(FT_y/(exp(2))));
            ax = min(idx_x); bx = max(idx_x);
            ay = min(idx_x); by = max(idx_y);
            MFD = max((x(bx)-x(ax)),(x(by)-x(ay))) ;
        case 'FWHM' 
            idx_x = find(FT_x > max(FT_x/2));
            idx_y = find(FT_y > max(FT_y/2));
            ax = min(idx_x); bx = max(idx_x);
            ay = min(idx_x); by = max(idx_y);
            MFD = max((x(bx)-x(ax)),(x(by)-x(ay))) ;
        case '4sigma'
            MFD = max(4*std(x,FT_x), 4*std(y,FT_y));
    end
end

            