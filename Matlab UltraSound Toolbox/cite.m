function cite

%CITE   How to cite MUST?
%   CITE suggests a list of references according to the functions that you
%   used for your project. Links to PDF articles are also provided in the
%   command window.
%
%
%   This function is part of MUST (Matlab UltraSound Toolbox).
%   MUST (c) 2020 Damien Garcia, LGPL-3.0-or-later
%
%   -- Damien Garcia -- 2022/04, last update 2024/01/18
%   website: <a
%   href="matlab:web('https://www.biomecardio.com')">www.BiomeCardio.com</a>


%-- Create a UI figure window
pos = get(groot,'defaultFigurePosition');
fig = uifigure;
fig.Position = [pos(1) pos(2) 310 280];
fig.Name = 'How to cite MUST?';

%-- Text
lbl = uilabel(fig);
lbl.Position = [20 230 300 30];
lbl.FontSize = 13;
lbl.Text = '<b>Please select what you used in MUST</b>:';
lbl.Interpreter = 'html';

%-- Create check boxes
%
%- DAS / DASMTX
cbx_DAS = uicheckbox(fig);
cbx_DAS.Position = [20 190 300 22];
cbx_DAS.Text = 'DAS, DASMTX, DAS3, or DASMTX3';
%
%- PFIELD
cbx_PFIELD = uicheckbox(fig);
cbx_PFIELD.Position = [20 160 200 22];
cbx_PFIELD.Text = 'PFIELD or PFIELD3';
%
%- SIMUS
cbx_SIMUS = uicheckbox(fig);
cbx_SIMUS.Position = [20 130 200 22];
cbx_SIMUS.Text = 'SIMUS or SIMUS3';
%
%- SPTRACK
cbx_SPTRACK = uicheckbox(fig);
cbx_SPTRACK.Position = [20 100 200 22];
cbx_SPTRACK.Text = 'SPTRACK';
%
%- Vector Doppler w/ PARAM.RXangle
cbx_RXangle = uicheckbox(fig);
cbx_RXangle.Position = [20 70 300 22];
cbx_RXangle.Text = 'The "RXangle" PARAM field for vector Doppler';
%
%- Automatic f-number
cbx_fnumber = uicheckbox(fig);
cbx_fnumber.Position = [20 40 300 22];
cbx_fnumber.Text = 'Automatic f-number: PARAM.fnumber = [ ]';

%-- DONE button
btn = uibutton(fig,'push');
btn.Text = 'Done';
btn.Position = [100 10 80 20];
btn.ButtonPushedFcn = @SelectionDone;

%--
function SelectionDone(~,~)
    RefList = {['D. Garcia, ',13,...
        '"Make the most of MUST, an open-source MATLAB UltraSound Toolbox," ',13,...
        '2021 IEEE International Ultrasonics Symposium (IUS), ',...
        '2021, pp. 1-4, ',13,...
        'doi: 10.1109/IUS52206.2021.9593605.']};
    urlList = {'https://www.biomecardio.com/publis/ius21.pdf'};
    
    %--- "Back to basics in ultrasound velocimetry"
    if cbx_SPTRACK.Value
        RefSPTRACK = ['V. Perrot and D. Garcia, ',13,...
            '"Back to basics in ultrasound velocimetry: tracking speckles by using a standard PIV algorithm," ',13,...
            '2018 IEEE International Ultrasonics Symposium (IUS), ',...'
            '2018, pp. 206-212, ',13,...
            'doi: 10.1109/ULTSYM.2018.8579665.'];
        RefList = [RefList,{''},RefSPTRACK];
        urlList = [urlList,'https://www.biomecardio.com/publis/ius18.pdf'];
    end
    
    %--- "So you think you can DAS?"
    if cbx_DAS.Value
        RefDAS = ['V. Perrot, M. Polichetti, F. Varray and D. Garcia, ',13,...
            '"So you think you can DAS? A viewpoint on delay-and-sum beamforming," ',13,...
            'Ultrasonics, 111, 2021, p. 106309,',13,...
            'doi: 10.1016/j.ultras.2020.106309.'];
        RefList = [RefList,{''},RefDAS];
        urlList = [urlList,'https://www.biomecardio.com/publis/ultrasonics21.pdf'];
    end
    %---

    %--- "SIMUS: An open-source simulator for medical ultrasound imaging"
    if cbx_PFIELD.Value || cbx_SIMUS.Value
        RefSIMUS1 = ['D. Garcia, ',13,...
            '"SIMUS: an open-source simulator for medical ultrasound imaging. Part I: theory & examples," ',13,...
            'Computer Methods and Programs in Biomedicine, 218, 2022, p. 106726, ',13,...
            'doi: 10.1016/j.cmpb.2022.106726.'];
        RefList = [RefList,{''},RefSIMUS1];
        urlList = [urlList,'https://www.biomecardio.com/publis/cmpb22.pdf'];
        RefSIMUS2 = ['A. Cigier, F. Varray and D. Garcia, ',13,...
            '"SIMUS: an open-source simulator for medical ultrasound imaging. Part II: comparison with four simulators," ',13,...
            'Computer Methods and Programs in Biomedicine, 218, 2022, p. 106774, ',13,...
            'doi: 10.1016/j.cmpb.2022.106774.'];
        RefList = [RefList,{''},RefSIMUS2];
        urlList = [urlList,'https://www.biomecardio.com/publis/cmpb22a.pdf'];
    end
    %---
    
    %--- "Color and vector flow imaging in parallel ultrasound"
    if cbx_RXangle.Value
        try
            %-- PMB 2018
            tmp = webread('https://scholar.google.com/citations?view_op=view_citation&citation_for_view=oS0jWX4AAAAJ:Cvh0bltMcLgC&hl=en');
            tmp = regexp(tmp,'Cited by ?(\d+)','tokens','once');
            cite1 = str2double(char(tmp));
            %-- TUFFC 2018
            tmp = webread('https://scholar.google.com/citations?view_op=view_citation&citation_for_view=oS0jWX4AAAAJ:It0W0vAlS5QC&hl=en');
            tmp = regexp(tmp,'Cited by ?(\d+)','tokens','once');
            cite2 = str2double(char(tmp));
        catch
            cite1 = rand(1);
            cite2 = rand(1);
        end
        if cite1>cite2
            RefRXangle = ['C. Madiena, J. Faurie, J. Porée and D. Garcia, ',13,...
                '"Color and vector flow imaging in parallel ultrasound with sub-Nyquist sampling," ',13,...
                'IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, ',...
                '65, 2018, pp. 795-802, ',13,...
                'doi: 10.1109/TUFFC.2018.2817885.'];
            urlList = [urlList,'https://www.biomecardio.com/publis/ieeeuffc18a.pdf'];
        else
            RefRXangle = ['S. Shahriari and D. Garcia, ',13,...
                '"Meshfree simulations of ultrasound vector flow imaging using smoothed particle hydrodynamics," ',13,...
                'Physics in Medicine and Biology, ',...
                '63, 2018, p. 205011, ',13,...
                'doi: 10.1088/1361-6560/aae3c3.'];
            urlList = [urlList,'https://www.biomecardio.com/publis/physmedbio18.pdf'];
        end
        RefList = [RefList,{''},RefRXangle];
        
    end
    %--- "Think twice before f-numbering"
    if cbx_fnumber.Value
        RefFNumber = ['D.Garcia and F. Varray, ',13,...
            '"Think twice before f-numbering," ',13,...
            'Ultrasonics, 2024, p. 107222, ',13,...
            'doi: 10.1016/j.ultras.2023.107222.'];
        RefList = [RefList,{''},RefFNumber];
        urlList = [urlList,'https://www.biomecardio.com/publis/ultrasonics24.pdf'];
    end
    %---    
    
    helpdlg(RefList,'Please cite:')
    
    RefList = RefList(~cellfun(@isempty,RefList));
    disp([13,'Thank you for making the most of MUST! Please cite:',13,...
        '------------------------------------------------------',13])
    for k = 1:numel(RefList)
        fprintf('%d) %s\n',k,RefList{k})
        fprintf('<a href = "%s">link to PDF</a>\n\n',urlList{k})
    end
    
    delete(fig)
end

end


