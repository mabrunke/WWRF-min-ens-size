"""
This is the code for Supplementary Section S1 Mean square skill score synthetic data
analysis in the following manuscript.

    Brunke et al. (2026). What is the minimum ensemble size that is sufficient to reduce
        uncertainty in assessing model skill in the West-WRF ensemble? Submitted to Monthly
        Weather Review.

The code is originally written by Leong Wai Siu.  If you have any questions, please send
an email to leongwaisiu@arizona.edu.

The code can be used to generate Figures S8 and S9.

"""
#!/usr/bin/env python
# Python library imports

# Third party imports
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# Local imports


class EnsembleMonteCarlo():

    def __init__(self, alpha, beta, nens=200, nsim=200, obs=0.5, mse_ref=0.5):
        self.alpha = alpha
        self.beta = beta
        self.nens = nens
        self.nsim = nsim
        self.obs = obs
        self.mse_ref = mse_ref

        # Initialize ensemble members
        ens = np.random.beta(self.alpha, self.beta, size=self.nens)
        ens_count, ens_bin = np.histogram(ens, bins=np.linspace(0,1,11))
        self.ens = ens
        self.ens_count = ens_count
        self.ens_bin = ens_bin

    def montecarlo(self):

        # Create working arrays
        ens_sample_mean = np.zeros((self.nsim, self.nens))
        mse_sim1 = np.zeros((self.nsim, self.nens))
        mse_sim2 = np.zeros((self.nsim, self.nens))
        msss_sim1 = np.zeros((self.nsim, self.nens))
        msss_sim2 = np.zeros((self.nsim, self.nens))

        # Default is 200 randomizations
        for j in np.arange(self.nsim):
            for i, nens_sample in enumerate(np.arange(self.nens)):
                # In each randomization, first time pick 1 member,
                # second time pick 2 members and so on.
                ens_sample_ind = np.random.choice(self.nens, nens_sample, replace=False)
                # Pick sample ensemble members
                ens_sample = self.ens[ens_sample_ind]
                ens_sample_mean[j,i] = np.mean(ens_sample)
                # Compute MSE using ensemble sample mean (Hans)
                mse_sim1[j,i] = np.square(ens_sample_mean[j,i] - self.obs)
                msss_sim1[j,i] = 1.0 - mse_sim1[j,i]/self.mse_ref
                # Compute MSE using all ensemble sample members (non-Hans)
                mse_sim2[j,i] = np.mean(np.square(ens_sample - self.obs))
                msss_sim2[j,i] = 1.0 - mse_sim2[j,i]/self.mse_ref

        # Mean of all randomized sampled mean; should approach ens_mean for all members
        ens_sim_mean = np.mean(ens_sample_mean, axis=0)
        # Mean of all randomized sampled MSE calculated from ensemble sample mean only
        mse_sim_mean1 = np.mean(mse_sim1, axis=0)
        # Mean of all randomized sampled MSE calculated from all samples
        mse_sim_mean2 = np.mean(mse_sim2, axis=0)
        # Mean of all randomized sampled MSSS calculated from ensemble sample mean only
        msss_sim_mean1 = np.mean(msss_sim1, axis=0)
        # Mean of all randomized sampled MSSS calculated from all samples
        msss_sim_mean2 = np.mean(msss_sim2, axis=0)

        return {
            "ens_sample_mean": ens_sample_mean,
            "mse_sim1": mse_sim1,
            "mse_sim2": mse_sim2,
            "msss_sim1": msss_sim1,
            "msss_sim2": msss_sim2,
            "ens_sim_mean": ens_sim_mean,
            "mse_sim_mean1": mse_sim_mean1,
            "mse_sim_mean2": mse_sim_mean2,
            "msss_sim_mean1": msss_sim_mean1,
            "msss_sim_mean2": msss_sim_mean2,
        }

def generate_ensemble_data(alpha, beta, nens):

    data = np.random.beta(alpha, beta, size=nens)

    return data

def plot_synthetic_msss_multiple_randomization(alpha, beta, nens=200, nrand=400, *,
                                               obs=0.5, mse_ref=0.5, fnum=1):
    """Generate synthetic data and plot different metrics."""
    ens = generate_ensemble_data(alpha, beta, nens)
    ens_count, ens_bin = np.histogram(ens, bins=np.linspace(0,1,11))
    # Get ensemble mean for all members
    ens_mean = np.mean(ens)

    # Create working arrays
    ens_sample_mean = np.zeros((nrand, nens))
    mse_rand1 = np.zeros((nrand, nens))
    mse_rand2 = np.zeros((nrand, nens))
    msss_rand1 = np.zeros((nrand, nens))
    msss_rand2 = np.zeros((nrand, nens))

    # Default is 400 randomizations
    for j in np.arange(nrand):
        # If 200 members, then do 200 trials
        for i, nens_sample in enumerate(np.arange(nens)):
            # In each randomization, first trial picks 1 member,
            # second trial picks 2 members and so on.
            ens_sample_ind = np.random.choice(nens, nens_sample, replace=False)
            # Pick sample ensemble members
            ens_sample = ens[ens_sample_ind]
            ens_sample_mean[j,i] = np.mean(ens_sample)
            # Compute MSE using ensemble sample mean (Hans)
            mse_rand1[j,i] = np.square(ens_sample_mean[j,i] - obs)
            msss_rand1[j,i] = 1.0 - mse_rand1[j,i]/mse_ref
            # Compute MSE using all ensemble sample members (non-Hans)
            mse_rand2[j,i] = np.mean(np.square(ens_sample - obs))
            msss_rand2[j,i] = 1.0 - mse_rand2[j,i]/mse_ref

    # Mean of all randomized sampled mean; should approach ens_mean for all members
    ens_rand_mean = np.mean(ens_sample_mean, axis=0)
    # Mean of all randomized sampled MSE calculated from ensemble sample mean only
    mse_rand_mean1 = np.mean(mse_rand1, axis=0)
    # Mean of all randomized sampled MSE calculated from all samples
    mse_rand_mean2 = np.mean(mse_rand2, axis=0)
    # Mean of all randomized sampled MSSS calculated from ensemble sample mean only
    # msss_rand_mean1 = np.mean(msss_rand1, axis=0)
    msss_rand_mean1 = 1.0 - mse_rand_mean1/mse_ref
    # Mean of all randomized sampled MSSS calculated from all samples
    # msss_rand_mean2 = np.mean(msss_rand2, axis=0)
    msss_rand_mean2 = 1.0 - mse_rand_mean2/mse_ref

    fig = plt.figure(figsize=[6.5,2.5], dpi=300, constrained_layout=False)
    gs = gridspec.GridSpec(1, 6, height_ratios=[1], width_ratios=[1,1,1,1,1,1],
        left=0.075, right=0.925, bottom=0.05, top=0.95, wspace=0.4, hspace=0.45)
    ax1, ax2, ax3, ax4, ax5, ax6 = [plt.subplot(gs[i]) for i in range(6)]
    axs = fig.axes

    # fig.suptitle(
    #     'Parameters: '
    #     + f'alpha = {alpha}; '
    #     + f'beta =  {beta}; '
    #     + f'nens =  {nens}; '
    #     + f'nrand =  {nrand}; '
    #     + f'obs =  {obs}; '
    #     + f'mse_ref =  {mse_ref}'
    #     )

    ens_count, ens_bin = np.histogram(ens, bins=np.linspace(0,1,11))
    ax1.hist(ens_bin[:-1],ens_bin,weights=ens_count/nens,
             lw=0.5,alpha=0.9,ec='k',fill=True,fc='none',rwidth=1.0)

    for j in np.arange(nrand):
        ax2.plot(np.arange(1,nens+1), ens_sample_mean[j,:], lw=0.6,c='lightgray')
        ax3.plot(np.arange(1,nens+1), mse_rand1[j,:], lw=0.6,c='lightgray')
        ax4.plot(np.arange(1,nens+1), mse_rand2[j,:], lw=0.6,c='lightgray')
        ax5.plot(np.arange(1,nens+1), msss_rand1[j,:], lw=0.6,c='lightgray')
        ax6.plot(np.arange(1,nens+1), msss_rand2[j,:], lw=0.6,c='lightgray')

    ax1.plot([ens_mean,ens_mean], [0,nens], lw=0.6,c='tab:red')
    ax1.plot([obs,obs], [0,nens], lw=0.6,linestyle='dashed', c='tab:red')
    ax2.plot([0,nens], [ens_mean,ens_mean], lw=0.6,c='tab:red')

    ax2.plot(np.arange(1,nens+1), ens_rand_mean, lw=0.6,c='k')
    ax3.plot(np.arange(1,nens+1), mse_rand_mean1, lw=0.6,c='k')
    ax4.plot(np.arange(1,nens+1), mse_rand_mean2, lw=0.6,c='k')
    ax5.plot(np.arange(1,nens+1), msss_rand_mean1, lw=0.6,c='k')
    ax6.plot(np.arange(1,nens+1), msss_rand_mean2, lw=0.6,c='k')

    plabels = ["(a)","(b)","(c)","(d)","(e)","(f)"]
    for i, ax in enumerate(axs):
        if i in [0]:
            ax.set_title('Ensemble PDF', pad=3)
            ax.set_xlabel('Ensemble value', labelpad=1)
            ax.set_ylabel('PDF', labelpad=1)
        elif i in [1]:
            ax.set_title('Ensemble mean', pad=3)
            ax.set_xlabel('Ensemble size', labelpad=1)
            ax.set_ylabel('Ensemble mean', labelpad=1)
        elif i in [2]:
            ax.set_title(r'MSE$_{ensm}$', pad=3)
            ax.set_xlabel('Ensemble size', labelpad=1)
            ax.set_ylabel('MSE', labelpad=1)
        elif i in [3]:
            ax.set_title(r'MSE$_{unc}$', pad=3)
            ax.set_xlabel('Ensemble size', labelpad=1)
            ax.set_ylabel('MSE', labelpad=1)
        elif i in [4]:
            ax.set_title(r'MSSS$_{ensm}$', pad=3)
            ax.set_xlabel('Ensemble size', labelpad=1)
            ax.set_ylabel('MSSS', labelpad=1)
        elif i in [5]:
            ax.set_title(r'MSSS$_{unc}$', pad=3)
            ax.set_xlabel('Ensemble size', labelpad=1)
            ax.set_ylabel('MSSS', labelpad=1)

        ax.text(0.15,0.94,plabels[i],
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes)
        ax.plot([0.05,0.05,0.25,0.25,0.05],
                [0.9,0.98,0.98,0.9,0.9],
                c="k",lw=0.4,
                transform=ax.transAxes)

    for i, ax in enumerate(axs):
        ax.tick_params(axis='both', which='major', length=3, pad=1, labelsize=5)
        ax.tick_params(axis='both', which='minor', length=2)
        # ax.xaxis.set_minor_locator(ticker.MultipleLocator(20))
        # ax.set_xlim(0.0,200.0)
        if i in [0]:
            ax.set_xticks(np.linspace(0.0,1.0,num=6))
            ax.set_ylim(0.0,0.5)
        elif i in [4,5]:
            ax.set_xticks(np.linspace(0.0,nens,num=6))
            ax.set_ylim(0.0,1.2)
        else:
            ax.set_xticks(np.linspace(0.0,nens,num=6))
            ax.set_ylim(0.0,1.2)
        ax.set_aspect(2.5/ax.get_data_ratio())

    outfile = f'/Volumes/craid/document/manuscript/journal/westwrf/figure/synthetic_msss_{fnum:02d}.png'
    plt.savefig(outfile)
    print(outfile)
    plt.clf()
    plt.close('all')
    plt.show()

def plot_synthetic_mse_single_randomization(alpha, beta, nens, *, fnum=1):

    ens = generate_ensemble_data(alpha, beta, nens)
    ens_mean = ens.mean()

    ens_count, ens_bin = np.histogram(ens, bins=np.linspace(0,1,11))
    obs = np.linspace(0,1,num=11)

    mse1 = [np.square(ens_mean - _) for _ in obs]
    mse2 = [np.square(ens - _).mean() for _ in obs]

    fig = plt.figure(figsize=[6.5,2.25], dpi=300, constrained_layout=False)
    gs = gridspec.GridSpec(1, 3, height_ratios=[1], width_ratios=[1,1,1],
        left=0.075, right=0.925, bottom=0.1, top=0.9, wspace=0.25, hspace=0.25)
    ax1, ax2, ax3 = [plt.subplot(gs[i]) for i in range(3)]
    axs = fig.axes

    print(ens_bin[:-1])
    print(ens_bin/nens)

    ax1.hist(ens_bin[:-1],ens_bin,weights=ens_count/nens,
             lw=0.5,alpha=0.9,ec='k',fill=True,fc='none',rwidth=1.0)
    ax2.plot(obs, mse1, lw=0.6,c='k')
    ax3.plot(obs, mse2, lw=0.6,c='k')
    # ax4.plot(obs, np.array(mse2)-np.array(mse1), lw=0.6,c='k')

    plabels = ["(a)","(b)","(c)"]
    for i, ax in enumerate(axs):
        if i in [0]:
            ax.set_title('Ensemble PDF', pad=2)
            ax.set_xlabel('Ensemble value', labelpad=1)
            ax.set_ylabel('PDF', labelpad=1)
        elif i in [1]:
            ax.set_title(r'MSE$_{ensm}$', pad=2)
            ax.set_xlabel('Obs', labelpad=1)
            ax.set_ylabel('MSE', labelpad=1)
        elif i in [2]:
            ax.set_title(r'MSE$_{unc}$', pad=2)
            ax.set_xlabel('Obs', labelpad=1)
            ax.set_ylabel('MSE', labelpad=1)
        # elif i in [3]:
        #     ax.set_title(r'(d) MSE$_2$ minus MSE$_1$', pad=2)
        #     ax.set_xlabel('Obs', labelpad=1)
        #     ax.set_ylabel('MSE difference', labelpad=1)

        ax.text(0.1,0.9,plabels[i],
                horizontalalignment='center',
                verticalalignment='center',
                transform=ax.transAxes)
        ax.plot([0.05,0.05,0.15,0.15,0.05],
                [0.85,0.95,0.95,0.85,0.85],
                c="k",lw=0.4,
                transform=ax.transAxes)

    for i, ax in enumerate(axs):
        ax.tick_params(axis='both', which='major', length=3, pad=1, labelsize=6)
        ax.tick_params(axis='both', which='minor', length=2)
        ax.set_xticks(np.linspace(0.0,1.0,num=6))
        # ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        # ax.set_xlim(0.0,1.0)
        if i in [0]:
            ax.set_ylim(0.0,0.5)
        elif i in [1,2]:
            ax.set_ylim(0.0,1.0)
        # elif i in [3]:
        #     ax.set_ylim(-0.2,0.2)
        ax.set_aspect(1.0/ax.get_data_ratio())

    outfile = f'/Volumes/craid/document/manuscript/journal/westwrf/figure/synthetic_mse_diff_{fnum:02d}.png'
    plt.savefig(outfile)
    print(outfile)
    plt.clf()
    plt.close('all')
    plt.show()

def main():

    alphas = [1, 3]
    betas = [4, 3]
    nens = 200
    for i, (alpha, beta) in enumerate(zip(alphas,betas)):
        plot_synthetic_mse_single_randomization(alpha, beta, nens, fnum=i+1)

    alpha_all = [1, 3]
    beta_all = [4, 3]
    # obs_all = [0.25,0.5,0.75]
    obs_all = [0.5,0.5]
    nens = 200
    nrand = 400
    mse_ref = 0.5
    # param_list = list(itertools.product(alpha_all,beta_all,obs_all))
    # for i, (alpha,beta,obs) in enumerate(param_list):
    for i, (alpha, beta,obs) in enumerate(zip(alpha_all,beta_all,obs_all)):
        plot_synthetic_msss_multiple_randomization(
            alpha, beta, nens, nrand,
            obs=obs, mse_ref=mse_ref, fnum=i+1
        )

if __name__ == "__main__":
    main()
