-- ec2k_axil_tb.vhd
-- The host program's view: drive ec2k_axil through its AXI4-Lite port,
-- replaying the WALK records of vectors_ecc2k130.txt across NENG engines.
-- Every record is loaded by register writes, every distinguished point is
-- read back through the queue and checked against the record, and at the
-- end the step and DP counters must agree with the records too.
--
-- With CDC=true the register block runs on its own clock of T_ENG_PS ps
-- behind the ec2k_axil_cdc bridge, as it does on F2 with the engine clock
-- above the shell's; the host side stays on the 10 ns clock.
--
--   ghdl -a --std=08 gf131_pkg.vhd gf131_tb_pkg.vhd gf2_kmul.vhd gf131_mul.vhd \
--        ec2k_batch_pipe.vhd ec2k_walker.vhd ec2k_axil.vhd ec2k_axil_cdc.vhd ec2k_axil_tb.vhd
--   ghdl -e --std=08 ec2k_axil_tb
--   ghdl -r --std=08 ec2k_axil_tb [-gNENG=1] [-gCDC=true -gT_ENG_PS=7000]

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;
use std.textio.all;

use work.gf131_pkg.all;
use work.gf131_tb_pkg.all;

entity ec2k_axil_tb is
  generic (
    VECTORS : string := "vectors_ecc2k130.txt";
    NENG    : natural := 2;
    ID_W    : natural := 4;
    LOG_W   : natural := 3;
    LOG_NB  : natural := 2;
    CDC     : boolean := false;
    T_ENG_PS : natural := 7000
  );
end entity;

architecture sim of ec2k_axil_tb is

  constant NWALK  : natural := NENG * 2 ** ID_W;
  constant MAXBLK : natural := 256;

  signal clk : std_logic := '0';
  signal rst : std_logic := '1';
  signal running : boolean := true;

  -- the register block's clock and reset: the host clock, or its own
  signal eclk : std_logic := '0';
  signal erst : std_logic := '1';

  signal awaddr, araddr, wdata, rdata : std_logic_vector(31 downto 0) := (others => '0');
  signal awvalid, awready, wvalid, wready, bvalid, bready : std_logic := '0';
  signal arvalid, arready, rvalid, rready : std_logic := '0';
  signal wstrb : std_logic_vector(3 downto 0) := "1111";
  signal bresp, rresp : std_logic_vector(1 downto 0);

  -- the register block's side of the bridge
  signal m_awaddr, m_araddr, m_wdata, m_rdata : std_logic_vector(31 downto 0);
  signal m_awvalid, m_awready, m_wvalid, m_wready, m_bvalid, m_bready : std_logic;
  signal m_arvalid, m_arready, m_rvalid, m_rready : std_logic;
  signal m_wstrb : std_logic_vector(3 downto 0);
  signal m_bresp, m_rresp : std_logic_vector(1 downto 0);

  type gf_array_t  is array (0 to MAXBLK - 1) of gf_t;
  type nat_array_t is array (0 to MAXBLK - 1) of natural;
  signal bx0, by0, bdx, bdy : gf_array_t := (others => (others => '0'));
  signal bk     : nat_array_t := (others => 0);
  signal nblk   : natural := 0;
  signal dpw    : natural := 0;
  signal loaded : boolean := false;

  signal errors : natural := 0;

  function word_of (v : gf_t; k : natural) return std_logic_vector is
    variable w : std_logic_vector(31 downto 0) := (others => '0');
  begin
    for b in 0 to 31 loop
      if 32 * k + b < M then
        w(b) := v(32 * k + b);
      end if;
    end loop;
    return w;
  end function;

begin

  clk <= not clk after 5 ns when running else '0';

  same_clock : if not CDC generate
    eclk <= clk;
    erst <= rst;
    m_awaddr <= awaddr; m_awvalid <= awvalid; awready <= m_awready;
    m_wdata <= wdata; m_wstrb <= wstrb; m_wvalid <= wvalid; wready <= m_wready;
    bresp <= m_bresp; bvalid <= m_bvalid; m_bready <= bready;
    m_araddr <= araddr; m_arvalid <= arvalid; arready <= m_arready;
    rdata <= m_rdata; rresp <= m_rresp; rvalid <= m_rvalid; m_rready <= rready;
  end generate;

  own_clock : if CDC generate
    eclk <= not eclk after (T_ENG_PS / 2) * 1 ps when running else '0';

    -- the engine-side reset releases a few of its clocks after the host's
    erst_sync : process (eclk)
      variable pipe : std_logic_vector(2 downto 0) := "111";
    begin
      if rising_edge(eclk) then
        pipe := pipe(1 downto 0) & rst;
        erst <= pipe(2);
      end if;
    end process;

    bridge : entity work.ec2k_axil_cdc
      generic map (SYNC_FF => 3)
      port map (
        s_clk => clk, s_rst => rst,
        s_awaddr => awaddr, s_awvalid => awvalid, s_awready => awready,
        s_wdata => wdata, s_wstrb => wstrb, s_wvalid => wvalid, s_wready => wready,
        s_bresp => bresp, s_bvalid => bvalid, s_bready => bready,
        s_araddr => araddr, s_arvalid => arvalid, s_arready => arready,
        s_rdata => rdata, s_rresp => rresp, s_rvalid => rvalid, s_rready => rready,
        m_clk => eclk, m_rst => erst,
        m_awaddr => m_awaddr, m_awvalid => m_awvalid, m_awready => m_awready,
        m_wdata => m_wdata, m_wstrb => m_wstrb, m_wvalid => m_wvalid, m_wready => m_wready,
        m_bresp => m_bresp, m_bvalid => m_bvalid, m_bready => m_bready,
        m_araddr => m_araddr, m_arvalid => m_arvalid, m_arready => m_arready,
        m_rdata => m_rdata, m_rresp => m_rresp, m_rvalid => m_rvalid, m_rready => m_rready);
  end generate;

  dut : entity work.ec2k_axil
    generic map (NENG => NENG, ID_W => ID_W, LOG_W => LOG_W, LOG_NB => LOG_NB,
                 FLUSH_CLK => 16, CNT_W => 32, DP_WEIGHT => 56, DP_FIFO_W => 3,
                 CLK_KHZ => 333333)
    port map (
      clk => eclk, rst => erst,
      s_awaddr => m_awaddr, s_awvalid => m_awvalid, s_awready => m_awready,
      s_wdata => m_wdata, s_wstrb => m_wstrb, s_wvalid => m_wvalid, s_wready => m_wready,
      s_bresp => m_bresp, s_bvalid => m_bvalid, s_bready => m_bready,
      s_araddr => m_araddr, s_arvalid => m_arvalid, s_arready => m_arready,
      s_rdata => m_rdata, s_rresp => m_rresp, s_rvalid => m_rvalid, s_rready => m_rready);

  loader : process
    file f        : text;
    variable ln   : line;
    variable st   : file_open_status;
    variable kind : string(1 to 5);
    variable h    : string(1 to HEXW);
    variable ch   : character;
    variable ok   : boolean;
    variable iv   : integer;
    variable n    : natural := 0;
  begin
    file_open(st, f, VECTORS, read_mode);
    assert st = open_ok report "cannot open " & VECTORS severity failure;
    while not endfile(f) and n < MAXBLK loop
      readline(f, ln);
      if ln'length >= 5 then
        read(ln, kind, ok);
        if ok and kind = "WALK " then
          read(ln, iv, ok); dpw <= iv;                read(ln, ch, ok);
          read(ln, h, ok);  bx0(n) <= hex_to_gf(h);  read(ln, ch, ok);
          read(ln, h, ok);  by0(n) <= hex_to_gf(h);
          read(ln, iv, ok); bk(n) <= iv;              read(ln, ch, ok);
          read(ln, h, ok);  bdx(n) <= hex_to_gf(h);  read(ln, ch, ok);
          read(ln, h, ok);  bdy(n) <= hex_to_gf(h);
          n := n + 1;
        end if;
      end if;
    end loop;
    file_close(f);
    nblk <= n;
    wait for 1 ns;
    loaded <= true;
    report "loaded " & integer'image(n) & " WALK records";
    wait;
  end process;

  host : process
    type blk_of_t is array (0 to NWALK - 1) of integer;
    variable cur      : blk_of_t := (others => -1);
    variable next_blk : natural := 0;
    variable v, st    : std_logic_vector(31 downto 0);
    variable w, g     : natural;
    variable px, py   : gf_t;
    variable reports  : natural := 0;
    variable polls    : natural := 0;
    variable want_steps : natural := 0;
    variable geom     : std_logic_vector(31 downto 0);

    -- one AXI-Lite write; address and data are offered together, the slave
    -- may take them on different clocks
    procedure wr (addr : natural; data : std_logic_vector(31 downto 0)) is
      variable aw_done, w_done : boolean := false;
    begin
      awaddr  <= std_logic_vector(to_unsigned(addr, 32));
      awvalid <= '1';
      wdata   <= data;
      wvalid  <= '1';
      bready  <= '1';
      while not (aw_done and w_done) loop
        wait until rising_edge(clk);
        if not aw_done and awready = '1' then
          aw_done := true;
          awvalid <= '0';
        end if;
        if not w_done and wready = '1' then
          w_done := true;
          wvalid <= '0';
        end if;
      end loop;
      while bvalid = '0' loop
        wait until rising_edge(clk);
      end loop;
      bready <= '0';
    end procedure;

    procedure rd (addr : natural; data : out std_logic_vector(31 downto 0)) is
    begin
      araddr  <= std_logic_vector(to_unsigned(addr, 32));
      arvalid <= '1';
      rready  <= '1';
      loop
        wait until rising_edge(clk);
        exit when arready = '1';
      end loop;
      arvalid <= '0';
      loop
        wait until rising_edge(clk);
        exit when rvalid = '1';
      end loop;
      data := rdata;
      rready <= '0';
    end procedure;

    procedure load (gid : natural; blk : natural) is
    begin
      wr(16#020#, std_logic_vector(to_unsigned(gid, 32)));
      for k in 0 to 4 loop
        wr(16#024# + 4 * k, word_of(bx0(blk), k));
        wr(16#038# + 4 * k, word_of(by0(blk), k));
      end loop;
      wr(16#04C#, x"00000001");
      loop
        rd(16#008#, st);
        exit when st(0) = '0';
      end loop;
      cur(gid) := blk;
    end procedure;

  begin
    wait until loaded;
    assert dpw = 56 report "vector cutoff " & integer'image(dpw)
      & " does not match the DP_WEIGHT the block was built with" severity failure;
    assert nblk >= NWALK report "need at least " & integer'image(NWALK) & " WALK records"
      severity failure;
    for i in 0 to nblk - 1 loop
      want_steps := want_steps + bk(i);
    end loop;

    wait until rising_edge(clk);
    wait until rising_edge(clk);
    rst <= '0';
    for i in 1 to 8 loop                      -- the engine-side reset lags
      wait until rising_edge(clk);
    end loop;

    rd(16#000#, v);
    assert v = x"2C130001" report "bad MAGIC " & to_hstring(v) severity failure;
    rd(16#050#, v);
    assert to_integer(unsigned(v)) = 333333 report "bad CLOCK " & to_hstring(v) severity error;
    rd(16#00C#, geom);
    assert to_integer(unsigned(geom(7 downto 0))) = ID_W
       and to_integer(unsigned(geom(31 downto 24))) = NENG
       and to_integer(unsigned(geom(23 downto 16))) = 56
      report "bad GEOM " & to_hstring(geom) severity failure;

    -- register readback of the load registers before anything runs
    wr(16#024#, x"DEADBEEF");
    rd(16#024#, v);
    assert v = x"DEADBEEF" report "LD_X0 readback" severity error;
    wr(16#034#, x"FFFFFFFF");
    rd(16#034#, v);
    assert v = x"00000007" report "LD_X4 should keep 3 bits, got " & to_hstring(v) severity error;

    wr(16#004#, x"00000003");                      -- RUN + CLEAR

    for i in 0 to NWALK - 1 loop
      load(i, i);
    end loop;
    next_blk := NWALK;

    while reports < nblk loop
      rd(16#008#, st);
      polls := polls + 1;
      assert polls < 200000 report "timeout after " & integer'image(reports) & " reports"
        severity failure;
      if st(1) = '1' then
        rd(16#080#, v);  g := to_integer(unsigned(v));
        rd(16#084#, v);  w := to_integer(unsigned(v));
        for k in 0 to 4 loop
          rd(16#090# + 4 * k, v);
          for b in 0 to 31 loop
            if 32 * k + b < M then px(32 * k + b) := v(b); end if;
          end loop;
          rd(16#0A4# + 4 * k, v);
          for b in 0 to 31 loop
            if 32 * k + b < M then py(32 * k + b) := v(b); end if;
          end loop;
        end loop;
        wr(16#0B8#, x"00000001");

        if g >= NWALK or cur(g) < 0 then
          report "report from unknown or idle walk " & integer'image(g) severity error;
          errors <= errors + 1;
        else
          if w /= bk(cur(g)) then
            report "walk " & integer'image(g) & " record " & integer'image(cur(g))
                   & ": reported after " & integer'image(w) & " steps, want "
                   & integer'image(bk(cur(g))) severity error;
            errors <= errors + 1;
          end if;
          if px /= bdx(cur(g)) or py /= bdy(cur(g)) then
            report "walk " & integer'image(g) & " record " & integer'image(cur(g))
                   & ": point mismatch, got (" & gf_to_hex(px) & ", " & gf_to_hex(py)
                   & ") want (" & gf_to_hex(bdx(cur(g))) & ", " & gf_to_hex(bdy(cur(g))) & ")"
              severity error;
            errors <= errors + 1;
          end if;
          cur(g) := -1;
          if next_blk < nblk then
            load(g, next_blk);
            next_blk := next_blk + 1;
          end if;
        end if;
        reports := reports + 1;
      end if;
    end loop;

    -- let the flush finish, then check the counters
    for i in 1 to 64 loop
      wait until rising_edge(clk);
    end loop;
    rd(16#010#, v);
    rd(16#014#, st);
    if to_integer(unsigned(v)) /= want_steps or st /= x"00000000" then
      report "STEPS = " & integer'image(to_integer(unsigned(v))) & ", want "
             & integer'image(want_steps) severity error;
      errors <= errors + 1;
    end if;
    rd(16#018#, v);
    if to_integer(unsigned(v)) /= nblk then
      report "DPS = " & integer'image(to_integer(unsigned(v))) & ", want " & integer'image(nblk)
        severity error;
      errors <= errors + 1;
    end if;
    rd(16#01C#, v);
    if v /= x"00000000" then
      report "DROPPED = " & integer'image(to_integer(unsigned(v))) severity error;
      errors <= errors + 1;
    end if;
    rd(16#008#, st);
    if st(1) = '1' or st(2) = '1' then
      report "STATUS after drain = " & to_hstring(st) severity error;
      errors <= errors + 1;
    end if;

    -- CLEAR zeroes the counters, RUN=0 parks the engines
    wr(16#004#, x"00000002");
    rd(16#010#, v);
    if v /= x"00000000" then
      report "STEPS after CLEAR = " & to_hstring(v) severity error;
      errors <= errors + 1;
    end if;
    rd(16#004#, v);
    if v /= x"00000000" then
      report "CTRL after CLEAR = " & to_hstring(v) severity error;
      errors <= errors + 1;
    end if;

    wait until rising_edge(clk);
    if errors = 0 then
      if CDC then
        report "ec2k_axil_tb: " & integer'image(reports) & " distinguished points, "
               & integer'image(want_steps) & " steps through " & integer'image(NENG)
               & " engine(s) of " & integer'image(2 ** ID_W) & " walks, "
               & integer'image(polls) & " status polls, engines on a "
               & integer'image(T_ENG_PS) & " ps clock behind the bridge";
      else
        report "ec2k_axil_tb: " & integer'image(reports) & " distinguished points, "
               & integer'image(want_steps) & " steps through " & integer'image(NENG)
               & " engine(s) of " & integer'image(2 ** ID_W) & " walks, "
               & integer'image(polls) & " status polls";
      end if;
      report "ec2k_axil_tb: PASS";
    else
      report "ec2k_axil_tb: FAIL -- " & integer'image(errors) & " errors" severity failure;
    end if;
    running <= false;
    wait;
  end process;

end architecture;
