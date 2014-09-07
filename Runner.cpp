#include "Runner.h"

#include <vector>

#include "MyStrategy.h"

using namespace model;
using namespace std;

int main(int argc, char* argv[]) {
    if (argc == 4) {
        Runner runner(argv[1], argv[2], argv[3]);
        runner.run();
    } else {
        Runner runner("127.0.0.1", "31001", "0000000000000000");
        runner.run();
    }
    
    return 0;
}

Runner::Runner(const char* host, const char* port, const char* token)
        : remoteProcessClient(host, atoi(port)), token(token) {
}

void Runner::run() {
    remoteProcessClient.writeTokenMessage(token);
    int teamSize = remoteProcessClient.readTeamSizeMessage();
    remoteProcessClient.writeProtocolVersionMessage();
    Game game = remoteProcessClient.readGameContextMessage();

    vector<Strategy*> strategies;

    for (int strategyIndex = 0; strategyIndex < teamSize; ++strategyIndex) {
        Strategy* strategy = new MyStrategy;
        strategies.push_back(strategy);
    }

    PlayerContext* playerContext;

    while ((playerContext = remoteProcessClient.readPlayerContextMessage()) != NULL) {
        vector<Hockeyist> playerHockeyists = playerContext->getHockeyists();
        if ((int) playerHockeyists.size() != teamSize) {
            break;
        }

        vector<Move> moves;

        for (int hockeyistIndex = 0; hockeyistIndex < teamSize; ++hockeyistIndex) {
            Hockeyist playerHockeyist = playerHockeyists[hockeyistIndex];

            Move move;
            strategies[playerHockeyist.getTeammateIndex()]
                    ->move(playerHockeyist, playerContext->getWorld(), game, move);
            moves.push_back(move);
        }

        remoteProcessClient.writeMovesMessage(moves);

        delete playerContext;
    }

    for (int strategyIndex = 0; strategyIndex < teamSize; ++strategyIndex) {
        delete strategies[strategyIndex];
    }
}
