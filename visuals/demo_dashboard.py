"""Interactive demo dashboard for EV aggregation results.

Run with:
    python3 demo_dashboard.py

This script loads the economic_savings.csv dataset and allows you to
experiment with different aggregator and alignment parameters. Use the
slider to change the maximum number of flex offers considered and press
"Show Demo" to update the plots.
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons, Button

DATA_FILE = Path(__file__).resolve().parent.parent / "data" / "economic_savings.csv"

AGGREGATOR_LABELS = {0: "Normal", 1: "TEC", 2: "Other"}
ALIGNMENT_LABELS = {0: "start", 1: "balance", 2: "price"}


def load_data() -> pd.DataFrame:
    """Load the economic savings dataset with user-friendly labels."""
    df = pd.read_csv(DATA_FILE)
    df["aggregator_label"] = df["aggregator_type"].map(AGGREGATOR_LABELS).fillna("Unknown")
    df["alignment_label"] = df["alignment"].map(ALIGNMENT_LABELS).fillna("Unknown")
    return df


def summarize_selection(df: pd.DataFrame, aggregator: int, alignment: int, max_flex: int):
    """Filter and summarize the dataset for the current control values."""
    filtered = df[
        (df["aggregator_type"] == aggregator)
        & (df["alignment"] == alignment)
        & (df["NrOfFlexOffers"] <= max_flex)
    ]

    if filtered.empty:
        return None, None, None

    average = filtered[["baseline_cost", "aggregated_cost", "savings", "scenario_time"]].mean()
    savings_by_flex = (
        filtered.groupby("NrOfFlexOffers")["savings"].mean().reset_index()
    )
    return average, savings_by_flex, filtered


def build_dashboard(df: pd.DataFrame):
    """Create the matplotlib figure with controls and initial plots."""
    fig, (ax_costs, ax_savings) = plt.subplots(2, 1, figsize=(10, 8))
    plt.subplots_adjust(left=0.1, bottom=0.2, right=0.75, hspace=0.45)

    # Controls
    aggregator_ax = fig.add_axes([0.78, 0.65, 0.18, 0.15])
    alignment_ax = fig.add_axes([0.78, 0.45, 0.18, 0.15])
    slider_ax = fig.add_axes([0.15, 0.1, 0.55, 0.03])
    button_ax = fig.add_axes([0.78, 0.25, 0.18, 0.07])

    aggregator_selector = RadioButtons(
        aggregator_ax,
        labels=[f"{k}: {v}" for k, v in AGGREGATOR_LABELS.items()],
        active=0,
    )
    alignment_selector = RadioButtons(
        alignment_ax,
        labels=[f"{k}: {v}" for k, v in ALIGNMENT_LABELS.items()],
        active=0,
    )

    max_offers = int(df["NrOfFlexOffers"].max())
    slider = Slider(slider_ax, "Max flex offers", 1, max_offers, valinit=max_offers, valstep=1)

    button = Button(button_ax, "Show demo", color="#4CAF50", hovercolor="#66bb6a")

    def refresh(_=None):
        aggregator_choice = int(aggregator_selector.value_selected.split(":", 1)[0])
        alignment_choice = int(alignment_selector.value_selected.split(":", 1)[0])
        max_flex_choice = int(slider.val)

        averages, savings_curve, filtered = summarize_selection(
            df, aggregator_choice, alignment_choice, max_flex_choice
        )

        ax_costs.clear()
        ax_savings.clear()

        title = (
            f"Aggregator: {AGGREGATOR_LABELS.get(aggregator_choice)} | "
            f"Alignment: {ALIGNMENT_LABELS.get(alignment_choice)} | "
            f"<= {max_flex_choice} flex offers"
        )
        fig.suptitle(title, fontsize=14, y=0.96)

        if averages is None:
            ax_costs.text(0.5, 0.5, "No data for this selection", ha="center", va="center")
            ax_savings.set_visible(False)
            fig.canvas.draw_idle()
            return

        ax_savings.set_visible(True)

        # Cost comparison bars
        ax_costs.bar(["Baseline cost", "Aggregated cost"], averages[["baseline_cost", "aggregated_cost"]])
        ax_costs.set_ylabel("Average cost")
        ax_costs.set_title("Cost comparison")

        # Annotate savings and scenario time
        savings_pct = (1 - averages["aggregated_cost"] / averages["baseline_cost"]) * 100
        summary = (
            f"Average savings: {averages['savings']:.2f}\n"
            f"Scenario time: {averages['scenario_time']:.3f}s\n"
            f"Efficiency gain: {savings_pct:.1f}%"
        )
        ax_costs.text(1.05, averages["aggregated_cost"], summary, va="bottom")

        # Savings curve
        ax_savings.plot(
            savings_curve["NrOfFlexOffers"],
            savings_curve["savings"],
            marker="o",
            color="#1f77b4",
        )
        ax_savings.set_xlabel("Number of flex offers")
        ax_savings.set_ylabel("Average savings")
        ax_savings.set_title("Savings as fleet size grows")
        ax_savings.grid(True, linestyle="--", alpha=0.4)

        fig.canvas.draw_idle()

    button.on_clicked(refresh)

    # Trigger first render
    refresh()
    plt.show()


def main():
    df = load_data()
    build_dashboard(df)


if __name__ == "__main__":
    main()
