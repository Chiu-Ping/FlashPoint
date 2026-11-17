#pragma once

class AnalysisConfig {
public:
    enum class Mode { Full, Light };

    static void SetMode(Mode mode) {
        ModeStorage() = mode;
    }

    static Mode GetMode() {
        return ModeStorage();
    }

    static bool FullHits() {
        return ModeStorage() == Mode::Full;
    }

    static bool ResponseH2() {
        return ModeStorage() == Mode::Full;
    }

private:
    static Mode& ModeStorage() {
        static Mode mode = Mode::Full;
        return mode;
    }
};
