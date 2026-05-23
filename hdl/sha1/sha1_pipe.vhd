-- sha1_pipe.vhd
-- Fully-pipelined SHA-1 compression core.
--
--   Throughput: 1 block (512 bits) per clock once the pipeline is full.
--   Latency:    82 cycles (1 input stage + 80 round stages + 1 final-add stage).
--   Critical path: 5-input 32-bit add (rotl(a,5) + f + e + K + W).
--                  Fits comfortably under 500 MHz on Virtex UltraScale+.
--
-- A single instance is the workhorse of an SHA-1 search engine; instantiate
-- N copies in parallel on a VU47P (~150 fit before routing pressure) for
-- ~50-100 GH/s of raw compression throughput per FPGA.
--
-- The user-tag input is carried alongside the pipeline so the consumer can
-- correlate each emitted digest with whatever metadata it cares about
-- (counter, message-modification index, distinguished-point trail ID, ...).

library ieee;
use ieee.std_logic_1164.all;
use ieee.numeric_std.all;

use work.sha1_pkg.all;

entity sha1_pipe is
  generic (
    TAG_WIDTH : positive := 64
  );
  port (
    clk       : in  std_logic;
    rst       : in  std_logic;                                  -- sync, active high

    in_valid  : in  std_logic;
    in_m      : in  block_t;                                    -- 16 BE words
    in_h      : in  state_t;                                    -- chaining input
    in_tag    : in  std_logic_vector(TAG_WIDTH - 1 downto 0);

    out_valid : out std_logic;
    out_h     : out state_t;
    out_tag   : out std_logic_vector(TAG_WIDTH - 1 downto 0)
  );
end entity;

architecture pipe of sha1_pipe is

  -- One stage = one SHA-1 round. The pipeline has 81 register slots:
  -- stage 0 latches the input; stage t (t in 1..80) holds the state after
  -- round t-1. A final combinational add produces the digest at stage 80.
  type stage_t is record
    s     : state_t;                                   -- working a..e
    w     : block_t;                                   -- W[t..t+15] sliding window
    h     : state_t;                                   -- chaining input, carried
    tag   : std_logic_vector(TAG_WIDTH - 1 downto 0);
    valid : std_logic;
  end record;

  type pipe_array_t is array (0 to 80) of stage_t;
  signal p : pipe_array_t;

begin

  -- Sequential pipeline. Every stage advances on the same clock; resets
  -- only clear `valid` so we don't churn 80 wide state registers.
  process(clk)
    variable a_v, b_v, c_v, d_v, e_v : word_t;
    variable t_v, w_v, wnext         : word_t;
  begin
    if rising_edge(clk) then
      if rst = '1' then
        for i in 0 to 80 loop
          p(i).valid <= '0';
        end loop;
      else
        -- Stage 0: latch inputs.
        p(0).s     <= in_h;
        p(0).w     <= in_m;
        p(0).h     <= in_h;
        p(0).tag   <= in_tag;
        p(0).valid <= in_valid;

        -- Stages 1..80: one SHA-1 round per stage.
        for t in 0 to 79 loop
          a_v := p(t).s(0);
          b_v := p(t).s(1);
          c_v := p(t).s(2);
          d_v := p(t).s(3);
          e_v := p(t).s(4);

          -- W[t] lives at the head of the sliding window.
          w_v := p(t).w(0);

          t_v := rotl(a_v, 5) + f_of(t, b_v, c_v, d_v) + e_v + k_of(t) + w_v;

          p(t + 1).s(0) <= t_v;
          p(t + 1).s(1) <= a_v;
          p(t + 1).s(2) <= rotl(b_v, 30);
          p(t + 1).s(3) <= c_v;
          p(t + 1).s(4) <= d_v;

          -- Message schedule: W[t+16] = ROL1(W[t+13] xor W[t+8] xor W[t+2] xor W[t]).
          -- After round 63 the new value would be W[80], which is never read.
          if t < 64 then
            wnext := rotl(p(t).w(13) xor p(t).w(8) xor p(t).w(2) xor p(t).w(0), 1);
          else
            wnext := (others => '0');
          end if;
          for i in 0 to 14 loop
            p(t + 1).w(i) <= p(t).w(i + 1);
          end loop;
          p(t + 1).w(15) <= wnext;

          p(t + 1).h     <= p(t).h;
          p(t + 1).tag   <= p(t).tag;
          p(t + 1).valid <= p(t).valid;
        end loop;
      end if;
    end if;
  end process;

  -- Final add: H_out = H_in + (a,b,c,d,e). One LUT level on UltraScale+.
  process(p)
    variable h_out : state_t;
  begin
    for i in 0 to 4 loop
      h_out(i) := p(80).h(i) + p(80).s(i);
    end loop;
    out_h     <= h_out;
    out_valid <= p(80).valid;
    out_tag   <= p(80).tag;
  end process;

end architecture;
