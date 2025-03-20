# python3
# Let each square piece be defined by the four colors of its four edges, in the format (up,left,down,right).

import sys
import math

sys.setrecursionlimit(10**6)

def parse_input():
    data = sys.stdin.read().splitlines()
    pieces = []
    for line in data:
        line = line.strip()
        if not line:
            continue
        # remove surrounding parentheses and then split on comma
        if line.startswith('(') and line.endswith(')'):
            line = line[1:-1]
        colors = line.split(',')
        # We assume there are exactly four colors: (up,left,down,right)
        pieces.append(tuple(colors))
    return pieces

def valid_piece(piece, r, c, n, board):
    # piece is a tuple: (up, left, down, right)
    # Check border conditions:
    if r == 0 and piece[0] != "black":
        return False
    if r == n-1 and piece[2] != "black":
        return False
    if c == 0 and piece[1] != "black":
        return False
    if c == n-1 and piece[3] != "black":
        return False
    # Check adjacent neighbors:
    # Check left neighbor: its right edge must match this piece's left edge.
    if c > 0:
        left_piece = board[r][c-1]
        if left_piece[3] != piece[1]:
            return False
    # Check top neighbor: its down edge must match this piece's up edge.
    if r > 0:
        top_piece = board[r-1][c]
        if top_piece[2] != piece[0]:
            return False
    return True

def backtrack(pos, n, pieces, used, board):
    if pos == n * n:
        return True
    r, c = divmod(pos, n)
    for i in range(len(pieces)):
        if not used[i]:
            piece = pieces[i]
            if valid_piece(piece, r, c, n, board):
                board[r][c] = piece
                used[i] = True
                if backtrack(pos + 1, n, pieces, used, board):
                    return True
                used[i] = False
                board[r][c] = None
    return False

def main():
    pieces = parse_input()
    total = len(pieces)
    # Determine n: n^2 pieces are given.
    n = int(math.sqrt(total))
    if n * n != total:
        sys.exit("Error: number of pieces is not a perfect square.")
    
    # Prepare board and a used flag list.
    board = [[None for _ in range(n)] for _ in range(n)]
    used = [False] * total

    if backtrack(0, n, pieces, used, board):
        # Print the board in the required format.
        # Each row is printed on one line, with pieces separated by semicolons,
        # and no spaces anywhere.
        for row in board:
            line = ';'.join("({},{},{},{})".format(piece[0], piece[1], piece[2], piece[3]) for piece in row)
            print(line)
    else:
        print("No valid placement found.")

if __name__ == '__main__':
    main()