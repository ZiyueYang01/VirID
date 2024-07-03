import matplotlib.pyplot as plt
import seaborn as sns
import math
import matplotlib.ticker as ticker

def format_xaxis(x, pos):
    return '{:.0f}k'.format(x/1000)


def clustr_pic(tsv, out_png):
    sns.set()
    sns.set_theme(style='ticks')
    df = tsv

    df.dropna(axis=0, how='any', inplace=True)
    groups = df.groupby(df.RdRP_super_group)
    clustr = df["RdRP_super_group"].unique()
    fig = plt.figure(figsize=(16,10)) 
    row = math.ceil(len(clustr)/3.0)
    columns = 3 if len(clustr) > 3 else len(clustr)

    for i in range(len(clustr)):
        plt.rcParams['figure.dpi'] = 75 
        group_df = groups.get_group(clustr[i])
        group_df.loc[:, "qlen"] = group_df["qlen"].astype(int)

        ax = fig.add_subplot(row, columns, i + 1)
        
        palette_color = {"noval":"#e76f51","known":"#219ebc"}

        bins = list(map(float, range(0, int(df['qlen'].max()+2000), 2000)))
        thistable = sns.histplot(data=group_df, x="qlen", hue="Virus_type",
                         bins=bins, shrink=0.8, multiple='dodge',
                         legend=True, palette=palette_color)
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                plt.text(p.get_x() + p.get_width() / 2.,
                        height + 0.05 * height,
                        f'{int(height):,}',
                        ha="center")

        thistable.set_xticks(range(0,int(df['qlen'].max()+2000),1000))

        sns.despine(top=True, right=True, left=False, bottom=False)

        ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_xaxis))
        ax.tick_params(axis='x', labelrotation=0)
        ax.set_title(clustr[i],pad=20)

        plt.subplots_adjust(wspace=0.3,hspace=0.7)
        
    plt.tight_layout() 
    plt.savefig(out_png, bbox_inches='tight')
