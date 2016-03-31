function PAR_YMD = dates_from_MD_to_YMD(rot, PAR_MD)
%% pre
insert_year = @(DATE,MD)  datestr( datenum(year(DATE),month(MD),day(MD)), 'yyyy-mm-dd' );
insert_year2= @(YEARS,MD) datestr( datenum(YEARS,     month(MD),day(MD)), 'yyyy-mm-dd' );
%% main
int         = cell(1);
int{1}      = insert_year( rot, PAR_MD );
% recognize the passage to a new year:
%   (when the the current date is smaller than the previous one ==> -1,
%    which happens when the year changes!)
jumps       = [1;sign(diff(datenum(int{1})))]==-1;
years       = year( rot ) + cumsum(jumps);
PAR_YMD     = cellstr(insert_year2( years, PAR_MD ));
%%
return