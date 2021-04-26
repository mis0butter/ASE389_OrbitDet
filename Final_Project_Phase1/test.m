% HW 5 
% Junette Hsin 

% finalproj_EOM_batch; 

%% EKF - all observations 

% Start figure 
ftitle = 'Radial-Intrack-Crosstrack'; 
figh = figure('name', ftitle); 
    hold on; grid on; 
    legend 

STATIONS = 0; 
% for STATIONS = 1:3 
for DATA = 0:2

    if STATIONS == 0    % use all station data 
            Yobs_STA = LEO_DATA_Apparent; 
        elseif STATIONS == 1 
            Yobs_STA = Yobs_KJL; 
        elseif STATIONS == 2
            Yobs_STA = Yobs_DGO;     
        else % STATIONS == 3
            Yobs_STA = Yobs_ACB;
        end 

        et_obs     = Yobs_STA(:,2) + et_t0; 
        XSTM_prev  = XSTM0_batch; 

        iter       = 0; 
        P_prev     = inv(Lambda_batch);

        X_EKF      = []; 
        t_X_EKF    = []; 
        Y_prefit   = []; 
        Y_postfit  = []; 

        Lambda_mat = []; 
        sigma3     = []; 

        for i = 1:length(et_obs)

            % keep track of iterations 
            iter = iter + 1; 
            sprintf('iter = %d', iter)

            % Propagate state 
            if     i == 1 && et_obs(1) == et_t0; t_prop = et_obs(i); 
            elseif i == 1;                       t_prop = [et_t0 : 60 : et_obs(1) ]; 
            else;                                t_prop = [et_obs(i-1) : 60 : et_obs(i)]; 
            end

            % EKF. All data, range, or range-only 
            if DATA == 0
                [t_XSTM, XSTM, Xstar, Y_pre, Y_post, P, Lambda] = fn.EKF(Yobs_STA, XSTM_prev, nX, ... 
                    epochs(1), t_prop, options, Amat_fn, Ht_fn, P_prev); 
            elseif DATA == 1
            % EKF for range only 
                [t_XSTM, XSTM, Xstar, Y_pre, Y_post, P, Lambda] = fn.EKF_r(Yobs_STA, XSTM_prev, nX, ... 
                epochs(1), t_prop, options, Amat_fn, Ht_r_fn, P_prev); 
            else % DATA == 2
            % EKF for range-rate only 
                [t_XSTM, XSTM, Xstar, Y_pre, Y_post, P, Lambda] = fn.EKF_rr(Yobs_STA, XSTM_prev, nX, ... 
                epochs(1), t_prop, options, Amat_fn, Ht_rr_fn, P_prev); 
            end

            % save states from current iteration 
            if i == 1 && et_obs(1) == et_t0; X_EKF = XSTM(1:nX)'; 
            else;                            X_EKF = [X_EKF; XSTM(:, 1:nX)]; 
            end 

            t_X_EKF   = [t_X_EKF; t_XSTM]; 
            Y_prefit  = [Y_prefit; Y_pre]; 
            Y_postfit = [Y_postfit; Y_post]; 

            % update measurement for next iteration 
            XSTM_prev = [Xstar; STM0]; 
            P_prev    = P;

            % innovations covariance 
            Lambda_mat = [Lambda_mat; Lambda]; 

            % 3-sigma STD 
            if DATA == 0; sigma3 = [sigma3; sqrt(Lambda(1,1))*3, sqrt(Lambda(2,2))*3];  
            else;         sigma3 = [sigma3; sqrt(Lambda(1,1))*3];  
            end

        end 

        % Propagate to last period of time 
        t_prop = [ et_obs(end) : 60 : et_t0 + 60*60*24 ]; 
        [t_XSTM, XSTM] = ode45(@(t, XSTM) fn.EOM_STM(t, XSTM, Amat_fn, nX), t_prop, XSTM_prev, options); 
        Xi   = XSTM(end,1:nX)'; 
        STMi = XSTM(end,nX+1:end); 
        STMi = reshape(STMi, [nX nX]); 

        % Time update + process noise 
        dt    = t_prop(end) - t_prop(1); 
        Q     = diag( (10000e-10)^2 * [1 1 1] ); 
        Gamma = [diag( dt^2/2 * [1 1 1] ); diag([dt dt dt])]; 
        P_noise = Gamma * Q * Gamma'; 
        Pi_bar = STMi * P_prev * STMi' + P_noise; 

        X_EKF_prop   = [X_EKF; XSTM(:,1:nX)]; 
        t_X_EKF_prop = [t_X_EKF; t_XSTM]; 


        %% Plot satellite position 

        ftitle = 'JahSat Orbit'; 
        figure('name', ftitle); 
            plot3(XSTM_ref0(:,1), XSTM_ref0(:,2), XSTM_ref0(:,3)); hold on; grid on; 
            plot3(XSTM_batch(:,1), XSTM_batch(:,2), XSTM_batch(:,3)); 
            plot3(X_EKF(:,1), X_EKF(:,2), X_EKF(:,3)); hold on; grid on; 
            plot3(X_EKF(1,1), X_EKF(1,2), X_EKF(1,3), 'o')
            xlabel('x (km)'); ylabel('y (km)'); zlabel('z (km)'); 
            legend('initial', 'batch', 'EKF', 'prop') 
            title(ftitle)


        %% Calculate residuals 

        Ypre_KJL = [];      Ypre_DGO = [];      Ypre_ACB = []; 

        Ypost_KJL = [];     Ypost_DGO = [];     Ypost_ACB = []; 

        sigma3_KJL = [];    sigma3_DGO = [];    sigma3_ACB = []; 

        % Extract states that correspond with station measurements 
        for i = 1:length(Y_postfit)

            ti  = Y_postfit(i, 2); 
            ti  = ti + et_t0; 
            i_X = find(t_X_EKF == ti); 

            i_STA = Y_postfit(i, 1); 

            if i_STA == 1
                Ypre_KJL   = [Ypre_KJL; Y_prefit(i,:)]; 
                Ypost_KJL  = [Ypost_KJL; Y_postfit(i,:)]; 
                sigma3_KJL = [sigma3_KJL; sigma3(i,:)]; 
            elseif i_STA == 2
                Ypre_DGO   = [Ypre_DGO; Y_prefit(i,:)]; 
                Ypost_DGO  = [Ypost_DGO; Y_postfit(i,:)]; 
                sigma3_DGO = [sigma3_DGO; sigma3(i,:)]; 
            else
                Ypre_ACB   = [Ypre_ACB; Y_prefit(i,:)]; 
                Ypost_ACB  = [Ypost_ACB; Y_postfit(i,:)]; 
                sigma3_ACB = [sigma3_ACB; sigma3(i,:)]; 
            end 

        end 

        t_KJL = Yobs_KJL(:,2); 
        t_DGO = Yobs_DGO(:,2); 
        t_ACB = Yobs_ACB(:,2); 

        if STATIONS == 0
            if DATA == 0
                res_all_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, 'Kwajalein Residuals'); 
                res_all_plot_STA(t_DGO, Yobs_DGO, Ypre_DGO, Ypost_DGO, sigma3_DGO, 'Diego-Garcia Residuals'); 
                res_all_plot_STA(t_ACB, Yobs_ACB, Yobs_ACB, Yobs_ACB, sigma3_ACB, 'Arecibo Residuals'); 
            elseif DATA == 1
                res_r_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, 'Kwajalein Residuals (Range Only)')
                res_r_plot_STA(t_DGO, Yobs_DGO, Ypre_DGO, Ypost_DGO, sigma3_DGO, 'Diego-Garcia Residuals (Range Only)')
                res_r_plot_STA(t_ACB, Yobs_ACB, Ypre_ACB, Ypost_ACB, sigma3_ACB, 'Arecibo Residuals (Range Only)')
            else
                res_rr_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, 'Kwajalein Residuals (Range-Rate Only)')
                res_rr_plot_STA(t_DGO, Yobs_DGO, Ypre_DGO, Ypost_DGO, sigma3_DGO, 'Diego-Garcia Residuals (Range-Rate Only)')
                res_rr_plot_STA(t_ACB, Yobs_ACB, Ypre_ACB, Ypost_ACB, sigma3_ACB, 'Arecibo Residuals (Range-Rate Only)')
            end
        elseif STATIONS == 1
            res_all_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, 'Kwajalein Residuals'); 
        elseif STATIONS == 2
            res_all_plot_STA(t_DGO, Yobs_DGO, Ypre_DGO, Ypost_DGO, sigma3_DGO, 'Diego-Garcia Residuals'); 
        else % STATIONS == 3
            res_all_plot_STA(t_ACB, Yobs_ACB, Yobs_ACB, Yobs_ACB, sigma3_ACB, 'Arecibo Residuals'); 
        end 


        %% radial-intrack-cross-track frame transformation for best estimate 

        % radial-intrack-cross-track frame transformation for best estimate 
        if STATIONS == 0  && DATA == 0 
            T_best = fn.ECItoRSW_T(Xi); 
        end 

        % transform all measurements 
        X_RSW = []; 
        for i = 1:length(X_EKF_prop) 
            X_RSW = [ X_RSW; [T_best * X_EKF_prop(i, 1:3)']' ]; 
        end 

        if STATIONS == 0
            if DATA == 0
                X_ALL_RSW = X_RSW; 
                save('X_ALL_RSW.mat') 
            elseif DATA == 1
                X_R_RSW = X_RSW; 
                save('X_R_RSW.mat') 
            else
                X_RR_RSW = X_RSW; 
                save('X_RR_RSW.mat') 
            end
        elseif STATIONS == 1
            X_KJL_RSW = X_RSW; 
            save('X_KJL_RSW.mat') 
        elseif STATIONS == 2
            X_DGO_RSW = X_RSW; 
            save('X_DGO_RSW.mat') 
        else
            X_ACB_RSW = X_RSW; 
            save('X_ACB_RSW.mat') 
        end 
            
        figure(figh); 
            plot3(X_RSW(:,1), X_RSW(:,2), X_RSW(:,3)); 

%     end
end


%% Subfunctions 

function res_all_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, ftitle)

    [dpre_err_KJL, dpre_rms_KJL, vpre_err_KJL, vpre_rms_KJL] = calc_res_all(Yobs_KJL, Ypre_KJL); 
    [dpost_err_KJL, dpost_rms_KJL, vpost_err_KJL, vpost_rms_KJL] = calc_res_all(Yobs_KJL, Ypost_KJL); 

    plot_res_all(ftitle, t_KJL, sigma3_KJL, dpre_err_KJL, dpre_rms_KJL, vpre_err_KJL, vpre_rms_KJL, ... 
        dpost_err_KJL, dpost_rms_KJL, vpost_err_KJL, vpost_rms_KJL)

end 

function res_r_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, ftitle)

    [dpre_err_KJL, dpre_rms_KJL]   = calc_res_r(Yobs_KJL, Ypre_KJL); 
    [dpost_err_KJL, dpost_rms_KJL] = calc_res_r(Yobs_KJL, Ypost_KJL); 

    plot_res_r(ftitle, t_KJL, sigma3_KJL, dpre_err_KJL, dpre_rms_KJL, dpost_err_KJL, dpost_rms_KJL)

end 

function res_rr_plot_STA(t_KJL, Yobs_KJL, Ypre_KJL, Ypost_KJL, sigma3_KJL, ftitle)

    [vpre_err_KJL, vpre_rms_KJL]   = calc_res_rr(Yobs_KJL, Ypre_KJL); 
    [vpost_err_KJL, vpost_rms_KJL] = calc_res_rr(Yobs_KJL, Ypost_KJL); 

    plot_res_r(ftitle, t_KJL, sigma3_KJL, vpre_err_KJL, vpre_rms_KJL, vpost_err_KJL, vpost_rms_KJL)

end 

function [d_err_STA, d_rms_STA, v_err_STA, v_rms_STA] = calc_res_all(Yobs_STA, Ycalc_STA)

% Calculate residuals 
d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
d_rms_STA = rms(d_err_STA); 
v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,4); 
v_rms_STA = rms(v_err_STA); 

end 

function [d_err_STA, d_rms_STA] = calc_res_r(Yobs_STA, Ycalc_STA)

% Calculate residuals 
d_err_STA = Yobs_STA(:,3) - Ycalc_STA(:,3); 
d_rms_STA = rms(d_err_STA); 

end 

function [v_err_STA, v_rms_STA] = calc_res_rr(Yobs_STA, Ycalc_STA)

% Calculate residuals 
v_err_STA = Yobs_STA(:,4) - Ycalc_STA(:,3); 
v_rms_STA = rms(v_err_STA); 

end 

function plot_res_all(ftitle, t_STA, sigma3_STA, dpre_err_STA, dpre_rms_STA, vpre_err_STA, vpre_rms_STA, ... 
    dpost_err_STA, dpost_rms_STA, vpost_err_STA, vpost_rms_STA)

figure('name', ftitle); 
    subplot(4,1,1) 
        scatter(t_STA, dpre_err_STA); hold on; grid on; 
        title({sprintf('PREFIT range residuals (km): mean = %.3g, RMS = %.3g', mean(dpre_err_STA), dpre_rms_STA)} ); 
        ylabel('km')  
    subplot(4,1,2) 
        scatter(t_STA, vpre_err_STA); hold on; grid on; 
        title({sprintf('PREFIT range-rate residuals (km/s): mean = %.3g, RMS = %.3g', mean(vpre_err_STA), vpre_rms_STA)} ); 
        xlabel('Time (s)') 
        ylabel('km/s') 
    subplot(4,1,3) 
        scatter(t_STA, dpost_err_STA); hold on; grid on; 
        plot(t_STA, sigma3_STA(:,1), 'r'); 
        plot(t_STA, -sigma3_STA(:,1), 'r'); 
        title({sprintf('POSTFIT range residuals (km): mean = %.3g, RMS = %.3g', mean(dpost_err_STA), dpost_rms_STA)} ); 
        ylabel('km')  
    subplot(4,1,4) 
        scatter(t_STA, vpost_err_STA); hold on; grid on; 
        plot(t_STA, sigma3_STA(:,2), 'r'); 
        plot(t_STA, -sigma3_STA(:,2), 'r'); 
        title({sprintf('POSTFIT range-rate residuals (km/s): mean = %.3g, RMS = %.3g', mean(vpost_err_STA), vpost_rms_STA)} ); 
        xlabel('Time (s)') 
        ylabel('km/s') 
    sgtitle(ftitle); 

end 

function plot_res_r(ftitle, t_STA, sigma3_STA, dpre_err_STA, dpre_rms_STA, dpost_err_STA, dpost_rms_STA)

figure('name', ftitle); 
    subplot(2,1,1) 
        scatter(t_STA, dpre_err_STA); hold on; grid on; 
        title({sprintf('PREFIT range residuals (km): mean = %.3g, RMS = %.3g', mean(dpre_err_STA), dpre_rms_STA)} ); 
        ylabel('km')  
    subplot(2,1,2) 
        scatter(t_STA, dpost_err_STA); hold on; grid on; 
        plot(t_STA, sigma3_STA(:,1), 'r'); 
        plot(t_STA, -sigma3_STA(:,1), 'r'); 
        title({sprintf('POSTFIT range residuals (km): mean = %.3g, RMS = %.3g', mean(dpost_err_STA), dpost_rms_STA)} ); 
        ylabel('km')  
    sgtitle(ftitle); 

end 

function plot_res_rr(ftitle, t_STA, sigma3_STA, vpre_err_STA, vpre_rms_STA, vpost_err_STA, vpost_rms_STA)

figure('name', ftitle); 
    subplot(2,1,1) 
        scatter(t_STA, vpre_err_STA); hold on; grid on; 
        title({sprintf('PREFIT range-rate residuals (km/s): mean = %.3g, RMS = %.3g', mean(vpre_err_STA), vpre_rms_STA)} ); 
        ylabel('km')  
    subplot(2,1,2) 
        scatter(t_STA, vpost_err_STA); hold on; grid on; 
        plot(t_STA, sigma3_STA(:,1), 'r'); 
        plot(t_STA, -sigma3_STA(:,1), 'r'); 
        title({sprintf('POSTFIT range-rate residuals (km/s): mean = %.3g, RMS = %.3g', mean(vpost_err_STA), vpost_rms_STA)} ); 
        ylabel('km')  
    sgtitle(ftitle); 

end 










