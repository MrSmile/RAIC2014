#pragma once

#ifndef _REMOTE_PROCESS_CLIENT_H_
#define _REMOTE_PROCESS_CLIENT_H_

#include <string>
#include <vector>

#include "csimplesocket/ActiveSocket.h"
#include "model/Game.h"
#include "model/Move.h"
#include "model/PlayerContext.h"
#include "model/World.h"

enum MessageType {
    UNKNOWN_MESSAGE,
    GAME_OVER,
    AUTHENTICATION_TOKEN,
    TEAM_SIZE,
    PROTOCOL_VERSION,
    GAME_CONTEXT,
    PLAYER_CONTEXT,
    MOVES_MESSAGE
};

class RemoteProcessClient {
private:
    CActiveSocket socket;
	bool cachedBoolFlag;
	bool cachedBoolValue;

    model::Game readGame();
    void writeGame(const model::Game& game);
    std::vector<model::Game> readGames();
    void writeGames(const std::vector<model::Game>& games);
    model::Hockeyist readHockeyist();
    void writeHockeyist(const model::Hockeyist& hockeyist);
    std::vector<model::Hockeyist> readHockeyists();
    void writeHockeyists(const std::vector<model::Hockeyist>& hockeyists);
    model::Move readMove();
    void writeMove(const model::Move& move);
    std::vector<model::Move> readMoves();
    void writeMoves(const std::vector<model::Move>& moves);
    model::Player readPlayer();
    void writePlayer(const model::Player& player);
    std::vector<model::Player> readPlayers();
    void writePlayers(const std::vector<model::Player>& players);
    model::PlayerContext readPlayerContext();
    void writePlayerContext(const model::PlayerContext& playerContext);
    std::vector<model::PlayerContext> readPlayerContexts();
    void writePlayerContexts(const std::vector<model::PlayerContext>& playerContexts);
    model::Puck readPuck();
    void writePuck(const model::Puck& puck);
    std::vector<model::Puck> readPucks();
    void writePucks(const std::vector<model::Puck>& pucks);
    model::World readWorld();
    void writeWorld(const model::World& world);
    std::vector<model::World> readWorlds();
    void writeWorlds(const std::vector<model::World>& worlds);

    static void ensureMessageType(MessageType actualType, MessageType expectedType);

    signed char readEnum();
    void writeEnum(signed char value);
    std::string readString();
    void writeString(const std::string& value);
    bool readBoolean();
    void writeBoolean(bool value);
    int readInt();
    void writeInt(int value);
    long long readLong();
    void writeLong(long long value);
    double readDouble();
    void writeDouble(double value);
    std::vector<signed char> readBytes(unsigned int byteCount);
    void writeBytes(const std::vector<signed char>& bytes);

    static bool isLittleEndianMachine();
public:
    RemoteProcessClient(std::string host, int port);

    void writeTokenMessage(const std::string& token);
    int readTeamSizeMessage();
    void writeProtocolVersionMessage();
    model::Game readGameContextMessage();
    model::PlayerContext* readPlayerContextMessage();
    void writeMovesMessage(const std::vector<model::Move>& move);

    void close();

    ~RemoteProcessClient();
};

#endif
