#include "RemoteProcessClient.h"

#include <algorithm>
#include <cstdlib>

using namespace model;
using namespace std;

const bool LITTLE_ENDIAN_BYTE_ORDER = true;
const int INTEGER_SIZE_BYTES = sizeof(int);
const int LONG_SIZE_BYTES = sizeof(long long);

RemoteProcessClient::RemoteProcessClient(string host, int port)
        : cachedBoolFlag(false), cachedBoolValue(false) {
    socket.Initialize();
    socket.DisableNagleAlgoritm();

    if (!socket.Open((uint8*) host.c_str(), (int16) port)) {
        exit(10001);
    }
}

void RemoteProcessClient::writeTokenMessage(const string& token) {
    writeEnum(AUTHENTICATION_TOKEN);
    writeString(token);
}

int RemoteProcessClient::readTeamSizeMessage() {
    ensureMessageType((MessageType) readEnum(), TEAM_SIZE);
    return readInt();
}

void RemoteProcessClient::writeProtocolVersionMessage() {
    writeEnum(PROTOCOL_VERSION);
    writeInt(1);
}

Game RemoteProcessClient::readGameContextMessage() {
    ensureMessageType((MessageType) readEnum(), GAME_CONTEXT);
    return readGame();
}

PlayerContext* RemoteProcessClient::readPlayerContextMessage() {
    MessageType messageType = (MessageType) readEnum();
    if (messageType == GAME_OVER) {
        return NULL;
    }

    ensureMessageType(messageType, PLAYER_CONTEXT);

    if (!readBoolean()) {
        return NULL;
    }

    cachedBoolFlag = true;
    cachedBoolValue = true;

    PlayerContext* playerContext = new PlayerContext;
    *playerContext = readPlayerContext();
    return playerContext;
}

void RemoteProcessClient::writeMovesMessage(const vector<Move>& moves) {
    writeEnum(MOVES_MESSAGE);
    writeMoves(moves);
}

void RemoteProcessClient::close() {
    socket.Close();
}

Game RemoteProcessClient::readGame() {
    if (!readBoolean()) {
        exit(20001);
    }

    long long randomSeed = readLong();
    int tickCount = readInt();
    double worldWidth = readDouble();
    double worldHeight = readDouble();
    double goalNetTop = readDouble();
    double goalNetWidth = readDouble();
    double goalNetHeight = readDouble();
    double rinkTop = readDouble();
    double rinkLeft = readDouble();
    double rinkBottom = readDouble();
    double rinkRight = readDouble();
    int afterGoalStateTickCount = readInt();
    int overtimeTickCount = readInt();
    int defaultActionCooldownTicks = readInt();
    int swingActionCooldownTicks = readInt();
    int cancelStrikeActionCooldownTicks = readInt();
    int actionCooldownTicksAfterLosingPuck = readInt();
    double stickLength = readDouble();
    double stickSector = readDouble();
    double passSector = readDouble();
    int hockeyistAttributeBaseValue = readInt();
    double minActionChance = readDouble();
    double maxActionChance = readDouble();
    double strikeAngleDeviation = readDouble();
    double passAngleDeviation = readDouble();
    double pickUpPuckBaseChance = readDouble();
    double takePuckAwayBaseChance = readDouble();
    int maxEffectiveSwingTicks = readInt();
    double strikePowerBaseFactor = readDouble();
    double strikePowerGrowthFactor = readDouble();
    double strikePuckBaseChance = readDouble();
    double knockdownChanceFactor = readDouble();
    double knockdownTicksFactor = readDouble();
    double maxSpeedToAllowSubstitute = readDouble();
    double substitutionAreaHeight = readDouble();
    double passPowerFactor = readDouble();
    double hockeyistMaxStamina = readDouble();
    double activeHockeyistStaminaGrowthPerTick = readDouble();
    double restingHockeyistStaminaGrowthPerTick = readDouble();
    double zeroStaminaHockeyistEffectivenessFactor = readDouble();
    double speedUpStaminaCostFactor = readDouble();
    double turnStaminaCostFactor = readDouble();
    double takePuckStaminaCost = readDouble();
    double swingStaminaCost = readDouble();
    double strikeStaminaBaseCost = readDouble();
    double strikeStaminaCostGrowthFactor = readDouble();
    double cancelStrikeStaminaCost = readDouble();
    double passStaminaCost = readDouble();
    double goalieMaxSpeed = readDouble();
    double hockeyistMaxSpeed = readDouble();
    double struckHockeyistInitialSpeedFactor = readDouble();
    double hockeyistSpeedUpFactor = readDouble();
    double hockeyistSpeedDownFactor = readDouble();
    double hockeyistTurnAngleFactor = readDouble();
    int versatileHockeyistStrength = readInt();
    int versatileHockeyistEndurance = readInt();
    int versatileHockeyistDexterity = readInt();
    int versatileHockeyistAgility = readInt();
    int forwardHockeyistStrength = readInt();
    int forwardHockeyistEndurance = readInt();
    int forwardHockeyistDexterity = readInt();
    int forwardHockeyistAgility = readInt();
    int defencemanHockeyistStrength = readInt();
    int defencemanHockeyistEndurance = readInt();
    int defencemanHockeyistDexterity = readInt();
    int defencemanHockeyistAgility = readInt();
    int minRandomHockeyistParameter = readInt();
    int maxRandomHockeyistParameter = readInt();
    double struckPuckInitialSpeedFactor = readDouble();
    double puckBindingRange = readDouble();

    return Game(randomSeed, tickCount, worldWidth, worldHeight, goalNetTop, goalNetWidth, goalNetHeight, rinkTop,
            rinkLeft, rinkBottom, rinkRight, afterGoalStateTickCount, overtimeTickCount, defaultActionCooldownTicks,
            swingActionCooldownTicks, cancelStrikeActionCooldownTicks, actionCooldownTicksAfterLosingPuck, stickLength,
            stickSector, passSector, hockeyistAttributeBaseValue, minActionChance, maxActionChance,
            strikeAngleDeviation, passAngleDeviation, pickUpPuckBaseChance, takePuckAwayBaseChance,
            maxEffectiveSwingTicks, strikePowerBaseFactor, strikePowerGrowthFactor, strikePuckBaseChance,
            knockdownChanceFactor, knockdownTicksFactor, maxSpeedToAllowSubstitute, substitutionAreaHeight,
            passPowerFactor, hockeyistMaxStamina, activeHockeyistStaminaGrowthPerTick,
            restingHockeyistStaminaGrowthPerTick, zeroStaminaHockeyistEffectivenessFactor, speedUpStaminaCostFactor,
            turnStaminaCostFactor, takePuckStaminaCost, swingStaminaCost, strikeStaminaBaseCost,
            strikeStaminaCostGrowthFactor, cancelStrikeStaminaCost, passStaminaCost, goalieMaxSpeed, hockeyistMaxSpeed,
            struckHockeyistInitialSpeedFactor, hockeyistSpeedUpFactor, hockeyistSpeedDownFactor,
            hockeyistTurnAngleFactor, versatileHockeyistStrength, versatileHockeyistEndurance,
            versatileHockeyistDexterity, versatileHockeyistAgility, forwardHockeyistStrength, forwardHockeyistEndurance,
            forwardHockeyistDexterity, forwardHockeyistAgility, defencemanHockeyistStrength,
            defencemanHockeyistEndurance, defencemanHockeyistDexterity, defencemanHockeyistAgility,
            minRandomHockeyistParameter, maxRandomHockeyistParameter, struckPuckInitialSpeedFactor, puckBindingRange);
}

void RemoteProcessClient::writeGame(const Game& game) {
    writeBoolean(true);

    writeLong(game.getRandomSeed());
    writeInt(game.getTickCount());
    writeDouble(game.getWorldWidth());
    writeDouble(game.getWorldHeight());
    writeDouble(game.getGoalNetTop());
    writeDouble(game.getGoalNetWidth());
    writeDouble(game.getGoalNetHeight());
    writeDouble(game.getRinkTop());
    writeDouble(game.getRinkLeft());
    writeDouble(game.getRinkBottom());
    writeDouble(game.getRinkRight());
    writeInt(game.getAfterGoalStateTickCount());
    writeInt(game.getOvertimeTickCount());
    writeInt(game.getDefaultActionCooldownTicks());
    writeInt(game.getSwingActionCooldownTicks());
    writeInt(game.getCancelStrikeActionCooldownTicks());
    writeInt(game.getActionCooldownTicksAfterLosingPuck());
    writeDouble(game.getStickLength());
    writeDouble(game.getStickSector());
    writeDouble(game.getPassSector());
    writeInt(game.getHockeyistAttributeBaseValue());
    writeDouble(game.getMinActionChance());
    writeDouble(game.getMaxActionChance());
    writeDouble(game.getStrikeAngleDeviation());
    writeDouble(game.getPassAngleDeviation());
    writeDouble(game.getPickUpPuckBaseChance());
    writeDouble(game.getTakePuckAwayBaseChance());
    writeInt(game.getMaxEffectiveSwingTicks());
    writeDouble(game.getStrikePowerBaseFactor());
    writeDouble(game.getStrikePowerGrowthFactor());
    writeDouble(game.getStrikePuckBaseChance());
    writeDouble(game.getKnockdownChanceFactor());
    writeDouble(game.getKnockdownTicksFactor());
    writeDouble(game.getMaxSpeedToAllowSubstitute());
    writeDouble(game.getSubstitutionAreaHeight());
    writeDouble(game.getPassPowerFactor());
    writeDouble(game.getHockeyistMaxStamina());
    writeDouble(game.getActiveHockeyistStaminaGrowthPerTick());
    writeDouble(game.getRestingHockeyistStaminaGrowthPerTick());
    writeDouble(game.getZeroStaminaHockeyistEffectivenessFactor());
    writeDouble(game.getSpeedUpStaminaCostFactor());
    writeDouble(game.getTurnStaminaCostFactor());
    writeDouble(game.getTakePuckStaminaCost());
    writeDouble(game.getSwingStaminaCost());
    writeDouble(game.getStrikeStaminaBaseCost());
    writeDouble(game.getStrikeStaminaCostGrowthFactor());
    writeDouble(game.getCancelStrikeStaminaCost());
    writeDouble(game.getPassStaminaCost());
    writeDouble(game.getGoalieMaxSpeed());
    writeDouble(game.getHockeyistMaxSpeed());
    writeDouble(game.getStruckHockeyistInitialSpeedFactor());
    writeDouble(game.getHockeyistSpeedUpFactor());
    writeDouble(game.getHockeyistSpeedDownFactor());
    writeDouble(game.getHockeyistTurnAngleFactor());
    writeInt(game.getVersatileHockeyistStrength());
    writeInt(game.getVersatileHockeyistEndurance());
    writeInt(game.getVersatileHockeyistDexterity());
    writeInt(game.getVersatileHockeyistAgility());
    writeInt(game.getForwardHockeyistStrength());
    writeInt(game.getForwardHockeyistEndurance());
    writeInt(game.getForwardHockeyistDexterity());
    writeInt(game.getForwardHockeyistAgility());
    writeInt(game.getDefencemanHockeyistStrength());
    writeInt(game.getDefencemanHockeyistEndurance());
    writeInt(game.getDefencemanHockeyistDexterity());
    writeInt(game.getDefencemanHockeyistAgility());
    writeInt(game.getMinRandomHockeyistParameter());
    writeInt(game.getMaxRandomHockeyistParameter());
    writeDouble(game.getStruckPuckInitialSpeedFactor());
    writeDouble(game.getPuckBindingRange());
}

vector<Game> RemoteProcessClient::readGames() {
    int gameCount = readInt();
    if (gameCount < 0) {
        exit(20002);
    }

    vector<Game> games;
    games.reserve(gameCount);

    for (int gameIndex = 0; gameIndex < gameCount; ++gameIndex) {
        games.push_back(readGame());
    }

    return games;
}

void RemoteProcessClient::writeGames(const vector<Game>& games) {
    int gameCount = games.size();
    writeInt(gameCount);

    for (int gameIndex = 0; gameIndex < gameCount; ++gameIndex) {
        writeGame(games[gameIndex]);
    }
}

Hockeyist RemoteProcessClient::readHockeyist() {
    if (!readBoolean()) {
        exit(20003);
    }

    long long id = readLong();
    long long playerId = readLong();
    int teammateIndex = readInt();
    double mass = readDouble();
    double radius = readDouble();
    double x = readDouble();
    double y = readDouble();
    double speedX = readDouble();
    double speedY = readDouble();
    double angle = readDouble();
    double angularSpeed = readDouble();
    bool teammate = readBoolean();
    HockeyistType type = (HockeyistType) readEnum();
    int strength = readInt();
    int endurance = readInt();
    int dexterity = readInt();
    int agility = readInt();
    double stamina = readDouble();
    HockeyistState state = (HockeyistState) readEnum();
    int originalPositionIndex = readInt();
    int remainingKnockdownTicks = readInt();
    int remainingCooldownTicks = readInt();
    int swingTicks = readInt();
    ActionType lastAction = (ActionType) readEnum();
    int lastActionTick = readBoolean() ? readInt() : -1;

    return Hockeyist(id, playerId, teammateIndex, mass, radius, x, y, speedX, speedY, angle, angularSpeed, teammate,
            type, strength, endurance, dexterity, agility, stamina, state, originalPositionIndex,
            remainingKnockdownTicks, remainingCooldownTicks, swingTicks, lastAction, lastActionTick);
}

void RemoteProcessClient::writeHockeyist(const Hockeyist& hockeyist) {
    writeBoolean(true);

    writeLong(hockeyist.getId());
    writeLong(hockeyist.getPlayerId());
    writeInt(hockeyist.getTeammateIndex());
    writeDouble(hockeyist.getMass());
    writeDouble(hockeyist.getRadius());
    writeDouble(hockeyist.getX());
    writeDouble(hockeyist.getY());
    writeDouble(hockeyist.getSpeedX());
    writeDouble(hockeyist.getSpeedY());
    writeDouble(hockeyist.getAngle());
    writeDouble(hockeyist.getAngularSpeed());
    writeBoolean(hockeyist.isTeammate());
    writeEnum((signed char) hockeyist.getType());
    writeInt(hockeyist.getStrength());
    writeInt(hockeyist.getEndurance());
    writeInt(hockeyist.getDexterity());
    writeInt(hockeyist.getAgility());
    writeDouble(hockeyist.getStamina());
    writeEnum((signed char) hockeyist.getState());
    writeInt(hockeyist.getOriginalPositionIndex());
    writeInt(hockeyist.getRemainingKnockdownTicks());
    writeInt(hockeyist.getRemainingCooldownTicks());
    writeInt(hockeyist.getSwingTicks());
    writeEnum((signed char) hockeyist.getLastAction());
    if (hockeyist.getLastActionTick() == -1) {
        writeBoolean(false);
    } else {
        writeBoolean(true);
        writeInt(hockeyist.getLastActionTick());
    }
}

vector<Hockeyist> RemoteProcessClient::readHockeyists() {
    int hockeyistCount = readInt();
    if (hockeyistCount < 0) {
        exit(20004);
    }

    vector<Hockeyist> hockeyists;
    hockeyists.reserve(hockeyistCount);

    for (int hockeyistIndex = 0; hockeyistIndex < hockeyistCount; ++hockeyistIndex) {
        hockeyists.push_back(readHockeyist());
    }

    return hockeyists;
}

void RemoteProcessClient::writeHockeyists(const vector<Hockeyist>& hockeyists) {
    int hockeyistCount = hockeyists.size();
    writeInt(hockeyistCount);

    for (int hockeyistIndex = 0; hockeyistIndex < hockeyistCount; ++hockeyistIndex) {
        writeHockeyist(hockeyists[hockeyistIndex]);
    }
}

Move RemoteProcessClient::readMove() {
    if (!readBoolean()) {
        exit(20005);
    }

    Move move;

    move.setSpeedUp(readDouble());
    move.setTurn(readDouble());
    move.setAction((ActionType) readEnum());
    if (move.getAction() == PASS) {
        move.setPassPower(readDouble());
        move.setPassAngle(readDouble());
    } else if (move.getAction() == SUBSTITUTE) {
        move.setTeammateIndex(readInt());
    }

    return move;
}

void RemoteProcessClient::writeMove(const Move& move) {
    writeBoolean(true);

    writeDouble(move.getSpeedUp());
    writeDouble(move.getTurn());
    writeEnum((signed char) move.getAction());
    if (move.getAction() == PASS) {
        writeDouble(move.getPassPower());
        writeDouble(move.getPassAngle());
    } else if (move.getAction() == SUBSTITUTE) {
        writeInt(move.getTeammateIndex());
    }
}

vector<Move> RemoteProcessClient::readMoves() {
    int moveCount = readInt();
    if (moveCount < 0) {
        exit(20006);
    }

    vector<Move> moves;
    moves.reserve(moveCount);

    for (int moveIndex = 0; moveIndex < moveCount; ++moveIndex) {
        moves.push_back(readMove());
    }

    return moves;
}

void RemoteProcessClient::writeMoves(const vector<Move>& moves) {
    int moveCount = moves.size();
    writeInt(moveCount);

    for (int moveIndex = 0; moveIndex < moveCount; ++moveIndex) {
        writeMove(moves[moveIndex]);
    }
}

Player RemoteProcessClient::readPlayer() {
    if (!readBoolean()) {
        exit(20007);
    }

    long long id = readLong();
    bool me = readBoolean();
    std::string name = readString();
    int goalCount = readInt();
    bool strategyCrashed = readBoolean();
    double netTop = readDouble();
    double netLeft = readDouble();
    double netBottom = readDouble();
    double netRight = readDouble();
    double netFront = readDouble();
    double netBack = readDouble();
    bool justScoredGoal = readBoolean();
    bool justMissedGoal = readBoolean();

    return Player(id, me, name, goalCount, strategyCrashed, netTop, netLeft, netBottom, netRight, netFront, netBack,
            justScoredGoal, justMissedGoal);
}

void RemoteProcessClient::writePlayer(const Player& player) {
    writeBoolean(true);

    writeLong(player.getId());
    writeBoolean(player.isMe());
    writeString(player.getName());
    writeInt(player.getGoalCount());
    writeBoolean(player.isStrategyCrashed());
    writeDouble(player.getNetTop());
    writeDouble(player.getNetLeft());
    writeDouble(player.getNetBottom());
    writeDouble(player.getNetRight());
    writeDouble(player.getNetFront());
    writeDouble(player.getNetBack());
    writeBoolean(player.isJustScoredGoal());
    writeBoolean(player.isJustMissedGoal());
}

vector<Player> RemoteProcessClient::readPlayers() {
    int playerCount = readInt();
    if (playerCount < 0) {
        exit(20008);
    }

    vector<Player> players;
    players.reserve(playerCount);

    for (int playerIndex = 0; playerIndex < playerCount; ++playerIndex) {
        players.push_back(readPlayer());
    }

    return players;
}

void RemoteProcessClient::writePlayers(const vector<Player>& players) {
    int playerCount = players.size();
    writeInt(playerCount);

    for (int playerIndex = 0; playerIndex < playerCount; ++playerIndex) {
        writePlayer(players[playerIndex]);
    }
}

PlayerContext RemoteProcessClient::readPlayerContext() {
    if (!readBoolean()) {
        exit(20009);
    }

    std::vector<Hockeyist> hockeyists = readHockeyists();
    World world = readWorld();

    return PlayerContext(hockeyists, world);
}

void RemoteProcessClient::writePlayerContext(const PlayerContext& playerContext) {
    writeBoolean(true);

    writeHockeyists(playerContext.getHockeyists());
    writeWorld(playerContext.getWorld());
}

vector<PlayerContext> RemoteProcessClient::readPlayerContexts() {
    int playerContextCount = readInt();
    if (playerContextCount < 0) {
        exit(20010);
    }

    vector<PlayerContext> playerContexts;
    playerContexts.reserve(playerContextCount);

    for (int playerContextIndex = 0; playerContextIndex < playerContextCount; ++playerContextIndex) {
        playerContexts.push_back(readPlayerContext());
    }

    return playerContexts;
}

void RemoteProcessClient::writePlayerContexts(const vector<PlayerContext>& playerContexts) {
    int playerContextCount = playerContexts.size();
    writeInt(playerContextCount);

    for (int playerContextIndex = 0; playerContextIndex < playerContextCount; ++playerContextIndex) {
        writePlayerContext(playerContexts[playerContextIndex]);
    }
}

Puck RemoteProcessClient::readPuck() {
    if (!readBoolean()) {
        exit(20011);
    }

    long long id = readLong();
    double mass = readDouble();
    double radius = readDouble();
    double x = readDouble();
    double y = readDouble();
    double speedX = readDouble();
    double speedY = readDouble();
    long long ownerHockeyistId = readLong();
    long long ownerPlayerId = readLong();

    return Puck(id, mass, radius, x, y, speedX, speedY, ownerHockeyistId, ownerPlayerId);
}

void RemoteProcessClient::writePuck(const Puck& puck) {
    writeBoolean(true);

    writeLong(puck.getId());
    writeDouble(puck.getMass());
    writeDouble(puck.getRadius());
    writeDouble(puck.getX());
    writeDouble(puck.getY());
    writeDouble(puck.getSpeedX());
    writeDouble(puck.getSpeedY());
    writeLong(puck.getOwnerHockeyistId());
    writeLong(puck.getOwnerPlayerId());
}

vector<Puck> RemoteProcessClient::readPucks() {
    int puckCount = readInt();
    if (puckCount < 0) {
        exit(20012);
    }

    vector<Puck> pucks;
    pucks.reserve(puckCount);

    for (int puckIndex = 0; puckIndex < puckCount; ++puckIndex) {
        pucks.push_back(readPuck());
    }

    return pucks;
}

void RemoteProcessClient::writePucks(const vector<Puck>& pucks) {
    int puckCount = pucks.size();
    writeInt(puckCount);

    for (int puckIndex = 0; puckIndex < puckCount; ++puckIndex) {
        writePuck(pucks[puckIndex]);
    }
}

World RemoteProcessClient::readWorld() {
    if (!readBoolean()) {
        exit(20013);
    }

    int tick = readInt();
    int tickCount = readInt();
    double width = readDouble();
    double height = readDouble();
    std::vector<Player> players = readPlayers();
    std::vector<Hockeyist> hockeyists = readHockeyists();
    Puck puck = readPuck();

    return World(tick, tickCount, width, height, players, hockeyists, puck);
}

void RemoteProcessClient::writeWorld(const World& world) {
    writeBoolean(true);

    writeInt(world.getTick());
    writeInt(world.getTickCount());
    writeDouble(world.getWidth());
    writeDouble(world.getHeight());
    writePlayers(world.getPlayers());
    writeHockeyists(world.getHockeyists());
    writePuck(world.getPuck());
}

vector<World> RemoteProcessClient::readWorlds() {
    int worldCount = readInt();
    if (worldCount < 0) {
        exit(20014);
    }

    vector<World> worlds;
    worlds.reserve(worldCount);

    for (int worldIndex = 0; worldIndex < worldCount; ++worldIndex) {
        worlds.push_back(readWorld());
    }

    return worlds;
}

void RemoteProcessClient::writeWorlds(const vector<World>& worlds) {
    int worldCount = worlds.size();
    writeInt(worldCount);

    for (int worldIndex = 0; worldIndex < worldCount; ++worldIndex) {
        writeWorld(worlds[worldIndex]);
    }
}

void RemoteProcessClient::ensureMessageType(MessageType actualType, MessageType expectedType) {
    if (actualType != expectedType) {
        exit(10011);
    }
}

signed char RemoteProcessClient::readEnum() {
    return this->readBytes(1)[0];
}

void RemoteProcessClient::writeEnum(signed char value) {
    vector<signed char> bytes;
    bytes.push_back(value);
    this->writeBytes(bytes);
}

string RemoteProcessClient::readString() {
    int length = this->readInt();
    if (length == -1) {
        exit(10014);
    }

    vector<signed char> bytes = this->readBytes(length);
    return string((char*) (&bytes[0]), length);
}

void RemoteProcessClient::writeString(const string& value) {
    vector<signed char> bytes(value.size());
    
    memcpy(&bytes[0], value.c_str(), value.size());

    this->writeInt(static_cast<int>(bytes.size()));
    this->writeBytes(bytes);
}

bool RemoteProcessClient::readBoolean() {
    if (cachedBoolFlag) {
        cachedBoolFlag = false;
        return cachedBoolValue;
    }
    return this->readBytes(1)[0] != 0;
}

void RemoteProcessClient::writeBoolean(bool value) {
    vector<signed char> bytes;
    bytes.push_back((signed char) (value ? 1 : 0));
    this->writeBytes(bytes);
}

int RemoteProcessClient::readInt() {
    vector<signed char> bytes = this->readBytes(INTEGER_SIZE_BYTES);

    if (this->isLittleEndianMachine() != LITTLE_ENDIAN_BYTE_ORDER) {
        reverse(&bytes[0], &bytes[INTEGER_SIZE_BYTES - 1]);
    }

    int value;

    memcpy(&value, &bytes[0], INTEGER_SIZE_BYTES);

    return value;
}

void RemoteProcessClient::writeInt(int value) {
    vector<signed char> bytes(INTEGER_SIZE_BYTES);

    memcpy(&bytes[0], &value, INTEGER_SIZE_BYTES);

    if (this->isLittleEndianMachine() != LITTLE_ENDIAN_BYTE_ORDER) {
        reverse(&bytes[0], &bytes[INTEGER_SIZE_BYTES - 1]);
    }

    this->writeBytes(bytes);
}

long long RemoteProcessClient::readLong() {
    vector<signed char> bytes = this->readBytes(LONG_SIZE_BYTES);

    if (this->isLittleEndianMachine() != LITTLE_ENDIAN_BYTE_ORDER) {
        reverse(&bytes[0], &bytes[LONG_SIZE_BYTES - 1]);
    }

    long long value;

    memcpy(&value, &bytes[0], LONG_SIZE_BYTES);

    return value;
}

void RemoteProcessClient::writeLong(long long value) {
    vector<signed char> bytes(LONG_SIZE_BYTES);

    memcpy(&bytes[0], &value, LONG_SIZE_BYTES);

    if (this->isLittleEndianMachine() != LITTLE_ENDIAN_BYTE_ORDER) {
        reverse(&bytes[0], &bytes[LONG_SIZE_BYTES - 1]);
    }

    this->writeBytes(bytes);
}

double RemoteProcessClient::readDouble() {
    long long value = this->readLong();
    return *((double*) &value);
}

void RemoteProcessClient::writeDouble(double value) {
    this->writeLong(*((long long*) &value));
}

vector<signed char> RemoteProcessClient::readBytes(unsigned int byteCount) {
    vector<signed char> bytes(byteCount);
    unsigned int offset = 0;
    int receivedByteCount;

    while (offset < byteCount && (receivedByteCount = socket.Receive(byteCount - offset)) > 0) {
        memcpy(&bytes[offset], socket.GetData(), receivedByteCount);
        offset += receivedByteCount;
    }

    if (offset != byteCount) {
        exit(10012);
    }

    return bytes;
}

void RemoteProcessClient::writeBytes(const vector<signed char>& bytes) {
    vector<signed char>::size_type byteCount = bytes.size();
    unsigned int offset = 0;
    int sentByteCount;

    while (offset < byteCount && (sentByteCount = socket.Send((uint8*) &bytes[offset], byteCount - offset)) > 0) {
        offset += sentByteCount;
    }

    if (offset != byteCount) {
        exit(10013);
    }
}

bool RemoteProcessClient::isLittleEndianMachine() {
    union {
        uint16 value;
        unsigned char bytes[2];
    } test = {0x0201};

    return test.bytes[0] == 1; 
}

RemoteProcessClient::~RemoteProcessClient() {
    this->close();
}
